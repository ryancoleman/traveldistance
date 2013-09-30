c include file for pymol style color maps: RGB range 0. - 1.
c--------------------------------------------
          real*4 color(3,3)
	  real*4 grey(3)
	  real*4 blue(3)
	  real*4 cyan(3)
	  real*4 green(3)
	  real*4 noir(3)
	  real*4 orange(3)
	  real*4 red(3)
	  real*4 violet(3)
	  real*4 white(3)
	  real*4 yellow(3)
	  real*4 yellow_green(3)
	  real*4 cyan_green(3)
	  real*4 cyan_blue(3)
	  real*4 violet_blue(3)
	  real*4 violet_red(3)
	  real*4 color1(3),color2(3),color3(3)
c--------------------------------------------
c B/W neutral
        data white / 1.0, 1.0, 1.0 /
        data grey / 0.4, 0.4, 0.4 /
        data noir / 0.1, 0.1, 0.1 /
c primary
        data red / 1.0, 0.0, 0.0 /
        data green / 0.0, 1.0, 0.0 /
        data blue / 0.0, 0.0, 1.0 /
c secondary
        data yellow / 1.0, 1.0, 0.0 /
        data cyan / 0.0, 1.0, 1.0 /
        data violet / 1.0, 0.0, 1.0 /
c tertiary
        data orange / 1.0, 0.5, 0.0 /
        data yellow_green / 0.5, 1.0, 0.0 /
        data cyan_green / 0.0, 1.0, 0.5 /
        data cyan_blue / 0.0, 0.5, 1.0 /
        data violet_red / 1.0, 0.0, 0.5 /
        data violet_blue / 0.5, 0.0, 1.0 /
c--------------------------------------------
