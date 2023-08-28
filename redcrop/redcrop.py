#!/usr/bin/env python3
"""
Do what the README file says.
"""

import sys, os, os.path as op, pathlib, re, functools, typing, argparse
from PIL import Image, ImageEnhance, ExifTags

mm_per_inch = 25.4

class CropLengths(typing.NamedTuple):
    left: float
    top: float
    right: float
    bottom: float

    @classmethod
    def from_string(CropLengths, s:str):
        try:
            return CropLengths(*[float(a.strip()) for a in s.split(",")])
        except ValueError:
            raise ValueError("Syntax for crop lenghts is "
                             "“left,top,right,bottom” in mm.")

    def __str__(self):
        return ",".join([str(int(a)) for a in self])

    def coordinates(self, imgsize):
        """
        Return a tuple suitable to pass to PIL.Image.crop().
        """
        width, height = imgsize
        return ( self.left,
                 self.top,
                 width - (self.left + self.right),
                 height - (self.top + self.bottom), )

resolution_re = re.compile(r"(\d*\.\d+)x(\d*\.\d+)")
class Resolution(typing.NamedTuple):
    x: float # dots per mm
    y: float # dots per mm

    @classmethod
    def from_string(Resulution, s:str):
        match = resolution_re.match(s)
        if match is None:
            try:
                r = int(s)
            except ValueError:
                raise ValueError("Resolution syntax is either two numers "
                                 "separated by a “x” as in “XxY” or a single "
                                 "number for equal resolition in both "
                                 "dimensions.")
            else:
                r = r / mm_per_inch
                return Resulution(r, r)
        else:
            x, y = match.groups()
            return Resolution(float(x) / mm_per_inch, float(y) / mm_per_inch)

    @classmethod
    def from_exif(Resolution, pil_image):
        exif = pil_image.getexif()

        x = exif.get(ExifTags.Base.XResolution, 0)
        y = exif.get(ExifTags.Base.YResolution, 0)

        if not x or not y:
            raise ValueError("Can’t determin resolution from EXIF data.")

        unit = exif[ExifTags.Base.ResolutionUnit]
        if unit == 2: # INCHES
            x = x / mm_per_inch
            y = y / mm_per_inch
        elif unit == 3: # CM
            x = x / 10.0
            y = y / 10.0
        else:
            raise ValueError(f"Unknown value for EXIF ResolutionUnit: {unit}.")

        return Resolution(x, y)

    def mm2px_x(self, mm):
        return mm*self.x

    def mm2px_y(self, mm):
        return mm*self.y

    @property
    def dpi_x(self):
        return self.x * mm_per_inch

    @property
    def dpi_y(self):
        return self.y * mm_per_inch

    def crop_in_px(self, crop_lengths):
        """
        Convert crop lengths specified in mm to pixels.
        """
        return CropLengths( self.mm2px_y(crop_lengths.left),
                            self.mm2px_x(crop_lengths.top),
                            self.mm2px_y(crop_lengths.right),
                            self.mm2px_x(crop_lengths.bottom) )

class CmdLineTool(object):
    def __init__(self):
        parser = self.argparse()
        self.args = parser.parse_args()

        if self.args.resolution is None:
            # Open first image file and get EXIF data.
            first = self.args.infilepaths[0]
            image = Image.open(first)
            self.resolution = Resolution.from_exif(image)
        else:
            self.resolution = self.args.resolution

        # Die Breite des Fadens in px.
        self.thread_width = self.resolution.mm2px_x(0.6)

        self.idx_tpl_re = re.compile(self.args.sort_index_re)

    @classmethod
    def argparse(CmdLineTool):
        """
        Return an ArgumentParser object
        """
        parser = argparse.ArgumentParser( prog=op.basename(sys.argv[0]),
                                          description=__doc__)

        parser.add_argument("--outdir", "-o", type=pathlib.Path,
                            default=pathlib.Path("."),
                            help="Output directory for page graphics")
        parser.add_argument("--text-width", "-w", type=float, required=True,
                            help="Width of the text block in mm")
        parser.add_argument("--inner-margin", "-i", type=float, required=True,
                            help="Distance of the text block from the red line")
        parser.add_argument("--safety", "-s", type=float, required=True,
                            default=4.0,
                            help="Space added to the text width left and "
                            "right for safety.")
        parser.add_argument("--margin-top", "-t", type=float, required=True,
                            help="Distance of the text block from top")
        parser.add_argument("--margin-bottom", "-b", type=float, required=True,
                            help="Distance of the text block from bottom")

        parser.add_argument("--resolution", "-r", type=Resolution.from_string,
                            default=None,
                            help="Input image resolution. If not specified, "
                            "attempt to read it the EXIF info of the first "
                            "input image.")

        parser.add_argument("--sort-index-re",
                            default="(\d+)(?:_(\d+))?\.[a-zA-Z+]$",
                            help="Regular expression containing groups "
                            "that each will be converted to integers "
                            "and then sorted as a tuple.")

        parser.add_argument("--pre-crop", type=CropLengths.from_string,
                            default=None,
                            help="Before processing, crop the image by the "
                            "specified length (comma-separated, mm, in order "
                            "left, top, right, bottom).")
        parser.add_argument("--pre-rotate", type=int, default=None,
                            help="Before processing, rotate the image by the "
                            "specified number of degrees.")

        parser.add_argument("-output-intermediates", "-O", action="store_true",
                            default=False,
                            help="Save intermediate stages of processing "
                            "to appropriately named files in the current "
                            "directory. Terminate after the first input file.")

        parser.add_argument("--contrast", "-C", type=int, default=4,
                            help="ImageEnhance.Contrast value")
        parser.add_argument("--brightness", "-B", type=int, default=4,
                            help="ImageEnhance.Brightness value")
        parser.add_argument("--scale-factor", "-S", type=float, default=1.5,
                            help="Greyscale to monochrome upscale factor")


        parser.add_argument("infilepaths", type=pathlib.Path, nargs="+",
                            help="Input image files")

        return parser



    def sortkey(self, name):
        match = self.idx_tpl_re.match(name)
        if match is None:
            raise OSError(f"Can’t create sort index for “{name}”.")
        return match.groups()

    def red_count(self, img):
        h = img.histogram()
        # Add up the counts for the reddest of the pixels.
        return (   sum(h[220:255]) # Red high
                   + sum(h[256:356]) # Green low
                   + sum(h[512:612]) # Blue low
                )

    def find_red_line(self, img):
        # In a range from the middle of the image
        # pick the one 50px band with the fewest non-yellow
        # pixels in it. (My brain hurts).
        area = self.resolution.mm2px_x(12)
        step = int(self.thread_width/2)
        width, height = img.size
        start = int((width/2)-(area/2))
        stop  = int(start + area)

        result = []
        for left in range(start, stop, step):
            #print(left, end=" - ")
            column = img.crop( (left, 0, int(left+step), height) )
            count = self.red_count(column)
            result.append( (int(left + (step/2)), count), )

            #print(" => ", count)

        left, least = max(result, key=lambda tpl: tpl[1])
        #print("result =", left, least)

        return left

    def main(self):
        def pageno_generator():
            p = 1
            while True:
                yield p
                p += 1
        pagenos = pageno_generator()


        stepnos = pageno_generator()
        def save_intermediate_step(name, img):
            if self.args.output_intermediates:
                stepno = next(stepnos)
                filename = f"{stepno:02d}_{name}.jpg"
                img.save(filename)

        sizes = {}
        for counter, infilepath in enumerate(self.args.infilepaths):
            img = Image.open(infilepath)

            if self.args.pre_rotate is not None:
                if self.resolution.x != self.resolution.y:
                    print("WARNING: Rotating the image might require "
                          "reversing x- and y-resolution. This has to be done"
                          "explicitly on the command line.", file=sys.stderr)

                img = img.rotate(int(self.args.pre_rotate), expand=True)
                save_intermediate_step(
                    f"pre_rotate_{self.args.pre_rotate}", img)

            if self.args.pre_crop is not None:
                crop = self.resolution.crop_in_px(self.args.pre_crop)
                img = img.crop(crop.coordinates(img.size))
                save_intermediate_step(f"pre_crop_{self.args.pre_crop}", img)

            width, height = img.size
            red_line = self.find_red_line(img)

            # Width of the text block on the page in pixels.
            text_width = self.resolution.mm2px_x(self.args.text_width)

            # Distance of the text block from the red line.
            inner_margin = self.resolution.mm2px_x(self.args.inner_margin)

            # Space added to the text width left and right
            # for safety.
            safety = self.resolution.mm2px_x(self.args.safety)

            # Distance of the text block from top and bottom.
            margin_top = self.resolution.mm2px_x(self.args.margin_top)
            margin_bottom = self.resolution.mm2px_x(self.args.margin_bottom)


            # img.crop ( (left, top, right, bottom,) )
            top = margin_top - safety
            bottom = height - margin_bottom + safety

            left_right = red_line - inner_margin
            left = img.crop( (left_right - text_width - 2*safety,
                              top,
                              left_right + safety,
                              bottom,) )

            right_left = red_line + inner_margin
            right = img.crop( (right_left - safety,
                               top,
                               right_left + text_width + 2*safety,
                               bottom,) )

            for img in ( left, right, ):
                # Convert to monochrome
                img = img.convert("L")

                # Contrast
                contrast = ImageEnhance.Contrast(img)
                img = contrast.enhance(self.args.contrast)

                #img.save("contrast.tif")

                # Brightness
                brightness = ImageEnhance.Brightness(img)
                img = brightness.enhance(self.args.brightness)

                #img.save("brightness.tif")

                scale_factor = self.args.scale_factor

                w, h = img.size
                img = img.resize( (int(w*scale_factor), int(h*scale_factor),) )

                # to Monochrome
                img = img.convert("1")

                # Save.
                fn = "%03d_%s.tif" % (next(pagenos), infilepath.stem,)
                outpath = pathlib.Path(self.args.outdir, fn)
                img.save(outpath,
                         dpi=(int(self.resolution.dpi_x*scale_factor),
                              int(self.resolution.dpi_y*scale_factor),))

            if self.args.output_intermediates:
                break

            print(counter, end="\r")
        print()
