# Finding chart for astronomical observations
This is a small python code to produce finding charts
STSCI archive. User can specify which band (survey) for the finding chart.
To run it, change the respective parameters in the init.yml file

Code allows to parse data from SIMBAD given a list of star names (from a list or a file) or
from a file containing the names, right ascension and declination of a list of stars.

The finding charts are obtained from STSCI archive. User can specify which band (survey) for the finding chart. 

Allows specification of the field of view.

**Examples for use**

0) `finding_chart.py -h` (or `-help`, `--h`, `--help`)
   - prints a small help screen
1) `finding_chart.py` 
   - reads the (standard) file `init.yml` and generate finding charts according to specification therein
2) `finding_chart.py init_starlist.yml` 
   - generates finding charts for the stars in the list specified in `init_starlist.yml` and outputs them in `./fc_my_starlist/`
3) `finding_chart.py init_starlist_file.yml` 
   - generates finding charts for the stars in the file `my_starlist.txt` and outputs them in `./fc_my_starlist_file/`
4) `finding_chart.py init_my_starnames_and_coords.yml`  
   - generates finding charts for the stars in the file `my_star_coords_txt` and outputs them in `./fc_my_star_coords_file/`
   
Use at your own peril and have fun observing

written by T. Maedler
