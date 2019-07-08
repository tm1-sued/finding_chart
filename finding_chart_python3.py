#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# written by Thomas MAEDLER
# 2019, July 6
#
# use at your own peril
#
import aplpy #for astro plotting pip install aplpy
import astropy.io.fits as pyfits
import matplotlib.pyplot  as plt
import numpy as np
import re
import os
from astroquery.simbad import Simbad
from astropy.io import ascii
from astropy.table import Table
from termcolor import colored
import astropy.io.fits as fits
import sys


def i_par():
    '''
    initialisation routine
    '''
    ip = {}
    
    #output folder for the finding charts
    ip['out_folder'] = './fc_my_starlist'
    #ip['out_folder'] = './fc_my_starlist_file'
    #ip['out_folder'] = './fc_my_star_coords_file'
    
    
    #type of input
    ip['fc_starname_input'] = 0 # provide star list according to simbad nomeclature in this initil file
    #ip['fc_starname_input'] = 1 # provide star list according to simbad nomeclature by a file containg a list of names
    #ip['fc_starname_input'] = 2 # provide user defined csv file with name and coordinates
    
    #parameter for fc_starname_input=0
    ip['star_list_for_simbad'] = ['HIP1','Cl* Blanco 1 WGL 28','Gaia DR2 5853498713160606720','Gliese 581']
    
    #parameter for fc_starname_input=1
    # path and name where is the list of stars
    ip['starlist_file'] = './my_starlist.txt' # single column file with the names of the stars retrieved from Simbad
    #name of output file for Simbad data (saved in out_folder)
    ip['ofile_simbad_data'] = 'my_starlist_Simbad_data.csv'
    
    # file name for option ip['fc_starname_input'] = 1 or 2
    ip['folder_stars_file'] = '.' # no / at the end
    
    
    
    ip['starfile_for_finding_chart'] = './my_star_coords_txt'
    
    # survey in which the finding chart is looked up
    # more info in
    # https://mast.stsci.edu/portal/Mashup/Clients/Mast/Portal.html
    # play with this not all catalogs give good finding charts
    ip['survey'] = 'dss1' # options all, dss1, poss2ukstu_red, poss2ukstu_ir, poss2ukstu_blue, poss1_blue, poss1_red, quickv, phase2_gsc2,phase2_gsc1
    
    ip['PI_name'] = 'G. Galilei'
    ip['Run_ID'] = 'Mike - 2019, July 6-8'
    
    
    
    #ip['star_file_head_name'] = 'source_id'
    ip['star_file_head_name'] = 'number_name'
    ip['star_file_head_ra'] = 'ra'
    ip['star_file_head_dec'] = 'dec'
    
    #ip['coord_format'] = 0 # decimal      (ra=23.1234 dec=-2.32456)
    #ip['coord_format'] = 1 # sexadecimal 1 (ra=23h01m01.23s dec=-02d32m04s56s) - SIMBAD format
    ip['coord_format'] = 2 # sexadecimal 2 (ra=23:01:01.23 dec=-02:32:04.56)
    
    #ip['diam_field_of_view_amin'] = 2.0 # this is the diameter in arcmin !! units right?
    ip['radius_field_of_view'] = 2.0 # this is the radius divided by 100 in deg 
    ip['radius_marking_circle_asec'] = 8.0
    
    
    return ip

#def read_file(path, file):
def read_file(path_file):
    '''
    reads a ascii table using the astropy frame work
    '''
    from astropy.io import ascii
    
    data = ascii.read(path_file)
    
    return data

def convert_ra_dec_to_decimal(ra_str, dec_str, form_choice):
    '''
    converts the string data for right ascension (ra_str) and declination (dec_str)
    that may be of the form
    
    form_choice = 0 : decimal       (e.g. ra=23.1234 dec=-2.32456)
    form_choice = 1 : sexadecimal 1 (e.g. ra=23h01m01.23s dec=-02d32m04s56s)
    form_choice = 2 : sexadecimal 1 (e.g. ra=23:01:01.23 dec=-02:32:04.56)
    
    to a decimals
    '''
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    
    ra_vec = []
    dec_vec = []
    ra_deg=float('NaN')
    dec_deg=float('NaN')
    
    if form_choice==0:
        return float(ra_str), float(dec_str)
    
    elif form_choice==1:
        if re.match('^(\d\d)h(\d\d)m(\d\d.\d+)s$',ra_str):
            ra_split = re.match('^(\d\d)h(\d\d)m(\d\d.\d+)s$',ra_str)
        else:
            ra_split = re.match('^(\d\d)h(\d\d)m(\d\d)s$',ra_str)
        
        if ra_split!=None:
            ra_vec = ra_split.groups()
        else:
            print('\n ra string does not comply with the format 00d00m00.000s')
            print('input was', repr(ra_str), '\n')
        
        if re.match('(-|\+)(\d\d)d(\d\d)m(\d\d.\d+)s$',dec_str):
            dec_split = re.match('(-|\+)(\d\d)d(\d\d)m(\d\d.\d+)s$',dec_str)
        else:
            dec_split = re.match('(-|\+)(\d\d)d(\d\d)m(\d\d)s$',dec_str)
        
        if dec_split!=None:
            dec_vec = dec_split.groups()
        else:
            print('\n dec string does not comply with the format 00d00m00.000s')
            print('input was', repr(dec_str), '\n')
    
    elif form_choice == 2:
        ra_split = re.match('^(\d\d):(\d\d):(\d\d.\d+)$',ra_str)
        if ra_split!=None:
            ra_vec = ra_split.groups()
        else:
            print('\n ra string does not comply with the format 00:00:00.000')
            print('input was', repr(ra_str), '\n')
        
        dec_split = re.match('(-|\+)(\d\d):(\d\d):(\d\d.\d+)$',dec_str)
        if dec_split!=None:
            dec_vec = dec_split.groups()
        else:
            print('\n dec string does not comply with the format 00:00:00.000')
            print('input was', repr(dec_str), '\n')
    else:
        print('RA/DEC format not specified, proviced format specifier', form_choice)
    
    if (len(ra_vec)!=0) and (len(dec_vec)!=0):
        m_ra = ':'.join(ra_vec)
        m_dec = dec_vec[0]+':'.join(dec_vec[1:])
        coords = SkyCoord(m_ra+' ' +m_dec, unit=(u.hourangle, u.deg))
        ra_deg = coords.ra.degree
        dec_deg = coords.dec.degree
    
    return ra_deg, dec_deg

def get_dss(ra, dec, starname, survey='dss1', radius = 2):
    '''
       uses the values for ra and dec to get a fits file of the sky from the 
       STSCI archive for
       survey \in {'all', 'dss1',...}
       in the window with height=width =radius with midpoint located at (ra,dec) 
       
    '''
    
    #  print "radius is in arcmins survey: 'all', 'poss2ukstu', 'dss1' 'dss2r' 'dss2b'..  "
    print ('getting fits file for '+ starname)
    print ('in survey bands: ', survey)
    
    url='http://archive.stsci.edu/cgi-bin/dss_search?v=%s&r=%fd&d=%fd&e=J2000&h=%f&w=%f&f=fits&c=none&fov=NONE&v3='%(survey,ra,dec,radius,radius)
    
    
    try:
        fitsdata = fits.open(url)
    except OSError:
        print('Fits file at url' + url)
        print( 'has '+colored('errors ', 'red')+'for star ' +  colored(starname, 'red'))
        print('... omitting production of finding chart')
        return None
    
    fitsfigure = aplpy.FITSFigure(fitsdata)
    print('data from URL',url)
    
    fitsfigure.show_grayscale(invert=True) #invert=true makes the image to be in negative format
    return fitsfigure

def get_list_from_file(folder_file, omit_hashed=True):
    '''
       reads a list form an ascii file in the folder with name file
       omits any line 
    '''
    lst=[]
    with open((folder_file))  as f:
        for line in f:
            if line[0]!='#' or not omit_hashed:
                lst.extend([line.strip('\n')])
    return lst

def get_simbad_data(starnames):
    """
        queries the coordinates, V mag, proper motion and spectral type from a list of stars
        outputs a table with Simbad ID, coordinates, Vmag, proper motion and spectral type
    """
    sim_data = {}
    
    mSimbad = Simbad()                    # only simbad ID & coord
    mSimbad.add_votable_fields('flux(V)') # add flux in V to table
    mSimbad.add_votable_fields('pm')      # add proper motion to table
    mSimbad.add_votable_fields('sp')      # add spectral type to table
    
    table_heads = ['my_starname', 'SimbadID', 'ra', 'dec', 'm_v', 'spectra_type', 'pm_ra', 'pm_dec']
    table_dtype = ['S50' for i in range(len(table_heads))]
    stable = Table(names = table_heads, dtype=table_dtype)
    
    for sn in starnames:
        print('...getting data for ', sn)
        newrow = [sn.replace(' ','_')]
        stardata = mSimbad.query_object(sn.replace('_', ' '))
        
        #print sn, stardata
        
        if stardata==None:
        	print(colored('could not get data for : %s' % sn, 'red'))
        else:
            print('newrow', newrow)
            #print('stardata', repr(stardata['MAIN_ID']))
            #print('stardata', stardata['MAIN_ID'][0].decode('UTF-8').split() )
            newrow.extend([' '.join(stardata['MAIN_ID'][0].decode('UTF-8').split()) ])
            
            ra = re.match('(\d\d) (\d\d) (\d\d\.\d+)$', stardata['RA'][0]).groups()
            try:
                dec = re.match('(-\d\d|\+\d\d) (\d\d) (\d\d\.\d+)$', stardata['DEC'][0]).groups()
            except:
                dec = re.match('(-\d\d|\+\d\d) (\d\d) (\d\d)$', stardata['DEC'][0]).groups()
                print (dec)
            
            m_ra  = ra[0]+u'h'+ra[1]+u'm'+ra[2]+u's'
            m_dec = dec[0]+u'd'+dec[1]+u'm'+dec[2]+u's'
            
            newrow.extend([m_ra, m_dec])
            newrow.extend([stardata['FLUX_V'][0]])
            newrow.extend([stardata['SP_TYPE'][0]])
            newrow.extend([stardata['PMRA'][0]])
            newrow.extend([stardata['PMDEC'][0]])
            
            stable.add_row(newrow)
    
    return stable

def label_fc(fc, sname, ra, dec, p):
    '''
        puts the labels around the finding chart
        adds an arrow to the plot pointing to the object
        adds the coordinates of the object in decmials
        adds a north/east directions
    '''
    if fc:
        #radius of field of view in arc sec
        #rad_fv_as = p['diam_field_of_view_amin']/60./2. # /60 to convert the arc min /2 to get radius
        rad_fv_as = p['radius_field_of_view']/60./2. # /60 to convert the arc min /2 to get radius
        #print('\nradius %f\n'% rad_fv_as)
        #print('\nradius.60 %f\n'% (rad_fv_as*60))
        fc.tick_labels.set_font(size='small')
        plt.tight_layout()
        fc.show_circles(ra,dec,radius=p['radius_marking_circle_asec']/3600.0,color='red')
        
        title_str = 'PI: %s ,   Run ID: %s,  Band: %s' % (p['PI_name'],\
        												  p['Run_ID'],\
        												  p['survey']) 
        
        plt.title(title_str ,fontsize=12)
        
        fc.tick_labels.set_xformat('dd.ddd') ## format to 2 decimals
        fc.tick_labels.set_yformat('dd.ddd') ## format to 2 decimals
        
        #this makes an orientation N/E coordinate system
        fc.add_label(ra-rad_fv_as*0.9,dec+rad_fv_as*0.95,'N',fontsize=10,color='blue')
        fc.add_label(ra-rad_fv_as*0.7,dec+rad_fv_as*0.75,'E',fontsize=10,color='blue')
        fc.show_arrows(ra-rad_fv_as*0.9,dec+rad_fv_as*0.75,dx=rad_fv_as*0.15,dy=0,width=0.2)
        fc.show_arrows(ra-rad_fv_as*0.9,dec+rad_fv_as*0.75,dx=0,dy=rad_fv_as*0.15,width=0.2)
        
        #this makes an arrow to the star
        fc.show_arrows(ra+rad_fv_as*0.45,dec-0.45*rad_fv_as, dx=-rad_fv_as*0.45 ,dy=rad_fv_as*0.45,width=0.2,color='red')
        
        #this displays the coordinates of the star
        label_str = r'$\alpha$ = %.6f   $\delta$= %.6f ' %(ra, dec)
        fc.add_label(ra+rad_fv_as*0.5, dec-0.5*rad_fv_as, sname,fontsize=15,color='red',weight='bold')
        fc.add_label(ra+rad_fv_as*0.4, dec-0.6*rad_fv_as,label_str,fontsize=15,color='red',weight='semibold')
        
        return fc

def make_fc_label_and_save(data, p):
    '''
        makes the finding chart, labels it and saves it under the specified location
        as a jpg
    '''
    
    if p['fc_starname_input'] in [0, 1]:
        coord_format = 1 # this is  how SIMBAD gives the coordinates
    else:
        coord_format = p['coord_format']
    
    for d in data:
        
        sname = str(d[p['star_file_head_name']])
        rastr = d[p['star_file_head_ra']]
        decstr =  d[p['star_file_head_dec']]
        
        ra_dec = convert_ra_dec_to_decimal(rastr, \
                                          decstr, \
                                          coord_format)
        ra = ra_dec[0]
        dec = ra_dec[1]
        
        im=get_dss(ra, dec, sname, survey=p['survey'], radius=p['radius_field_of_view'])
        #print(sname, ra_dec)
        
        fc  = label_fc(im, sname, ra, dec, p)
        if fc:
            im.save(p['out_folder']+'/'+sname.replace(' ', '_')+'.jpeg')


def main():
    '''
        main file to 
        a) parse the coordinates from SIMBAD using either a list 
           of star names (in SIMBAD nomenclature) like star_list = ['Rigel', 'HIP1']
           or from a file containg such list of star names 
        b) read an astropy table with colums star_name, right ascension, declination
        
        make the finding chart
    '''
    ip = i_par()
    
    if not os.path.exists(ip['out_folder']):
        os.makedirs(ip['out_folder'])
    
    
    
    if ip['fc_starname_input'] in [0, 1]:
        
        if ip['fc_starname_input']==0:
            
            star_list = ip['star_list_for_simbad']
        else:
            
            star_list = get_list_from_file(ip['starlist_file'])
        
        sim_data = get_simbad_data(star_list)
        ascii.write(sim_data, os.path.join(ip['out_folder'], ip['ofile_simbad_data']), format='csv')
        
        print('...Simbad data written to ', os.path.join(ip['out_folder'], ip['ofile_simbad_data']))
        # reset labels 
        ip['star_file_head_name'] = 'my_starname'
        ip['star_file_head_ra'] = 'ra'
        ip['star_file_head_dec'] = 'dec'
        
        make_fc_label_and_save(sim_data, ip)
        
    elif ip['fc_starname_input']==2:
        print('reading from my file')
        #data = read_file(ip['folder_stars_file'], ip['starfile_for_finding_chart'])
        data = read_file(ip['starfile_for_finding_chart'])
        
        print('keys of data file', data.keys())
        print('here')
        make_fc_label_and_save(data, ip)
        
    else:
        
        print("option 'fc_starname_input'=", ip['fc_starname_input'],)
        print('not available')
        print('...exit')

if __name__ == "__main__":
    main()