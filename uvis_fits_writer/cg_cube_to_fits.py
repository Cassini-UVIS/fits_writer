'''
Created on 29 Nov 2021

@author: Josh Elliott : joshua.elliott@lasp.colorado.edu

TODO: Old cubes vs new cubes.  Need to handle both.  For example, the cal factor matrix.
TODO: Handle fits, sav and PDS files
TODO: Handle cube files with only one datastruct (hdu)
TODO: Discussion on why original cubes removed ymin offset from all images and backplanes?  
      Shouldn't we maintain pixel position throughout?  For now, I'm just copying
      things over, leaving the images at y=0
    
'''

import numpy as np
from numpy import genfromtxt
from astropy.io import fits
from astropy.table import Table
from uvis_template_factory import UVISFitsTemplateFactory
from pathlib import Path
from datetime import datetime
import pvl
from spiceypy import spiceypy as cspice
import subprocess
from planetarypy.spice.kernels import list_kernels_for_day

# TODO: Move this to a general spice function library.
def get_orbit_number(sclk_time):
    spice_dir = Path('..') / 'kernels'
    orb_file = spice_dir / 'orb' / 'cas_v40.orb'
    orb_numbers = []
    orb_sclk_times = []
    with open(orb_file) as file:
        lines = file.readlines()[2:]
        for line in lines:
            parts = [p for p in line.split(' ') if p != '']
            orb_numbers.append(np.int(parts[0]))
            orb_sclk_times.append(np.double(parts[5].split('/')[1]))
            
    orb_numbers = np.array(orb_numbers)
    orb_sclk_times = np.array(orb_sclk_times)
    
    w = np.where(sclk_time >= orb_sclk_times)
    return orb_numbers[w[0][-1]]

def cassini_uvis_euv_wavelengths(xbin):
    '''
    This function was originally written in IDL by the Cassini UVIS processing team. 
    The IDL source code can be found at https://github.com/Cassini-UVIS/tools
    Ported to Python by Emma Lieb : https://github.com/emmalieb/SALSA
    
    :param xbin:   Number of wavelength elements to bin together.
    :type xbin:    int
    
    :return:       Array of wavelength values of size 1024/xbin
    :rtype:        numpy.array, float64
    
    TODO: This should be moved to a general library of UVIS functions.
    '''
    
    RAD=180.0/np.pi
    D=1.E7/1371.
    ALP=8.03/RAD+.00032-.0000565 
    BET=(np.linspace(0,1023,1024)-511.5)*0.025*.9987/300.0
    BET=np.arctan(BET)-1.19/RAD+.00032-.0000565
    LAM=D*(np.sin(ALP)+np.sin(BET))
    e_wavelength=LAM
    if xbin == 1:
        return e_wavelength
    
    e_wavelength=np.zeros(shape=(1024//xbin,))
    for k in range(0,1024//xbin):
        e_wavelength[k]=np.sum(LAM[k*xbin:(k+1)*xbin-1])/xbin
    
    return e_wavelength

def cassini_uvis_fuv_wavelengths(xbin):
    '''
    This function was originally written in IDL by the Cassini UVIS processing team. 
    The IDL source code can be found at https://github.com/Cassini-UVIS/tools
    Ported to Python by Emma Lieb : https://github.com/emmalieb/SALSA
    
    :param xbin:   Number of wavelength elements to bin together.
    :type xbin:    int
    
    :return:       Array of wavelength values of size 1024/xbin
    :rtype:        numpy.array, float64
    
    TODO: This should be moved to a general library of UVIS functions.
    '''
    RAD=180./np.pi
    D=1.E7/1066
    ALP=(9.22+.032)/RAD
    ALP=ALP+3.46465e-5
    BET=(np.linspace(0,1023,1024)-511.5)*0.025*0.99815/300.0
    BET=np.arctan(BET)+0.032/RAD+3.46465e-5
    lam=D*(np.sin(ALP)+np.sin(BET))
    e_wavelength=lam
    if xbin == 1:
        return e_wavelength
    
    e_wavelength=np.zeros(shape=(1024//xbin,))
    for k in range(0,1024//xbin):
        e_wavelength[k]=np.sum(lam[k*xbin:(k+1)*xbin])/xbin
    
    return e_wavelength

class CGCubeToFITS(object):
    '''
    This class maps an existing Cassini UVIS Cube file into the new FITS format
    given a template.
    '''
    
    def __init__(self):
        '''
        Constructor.
        '''
        pass
    
    def convert_to_fits(self, template_file, cube_file, output_file, pds_label_file, 
                        kernel_list=None):
        '''
        Convert a UVIS Cube file into the new FITS format.
        
        :param template_file:    Path to spreadsheet containing the template.
        :type template_file:     pathLib.Path object.
        :param cube_file:        Path to original cube file.
        :type cube_file:         pathLib.Path object.
        :param output_file:      Path to output FITS file.
        :type output_file:       pathLib.Path object.
        :param pds_label_file:   Path to PDS Label file for this observation
        :type pds_label_file:    pathLib.Path object.
        :param kernel_list:      List of kernel files used (for metadata purposes).
        :type kernel_list:       List of strings.
        '''
        
        # Create the template factory
        factory = UVISFitsTemplateFactory()
        self.factory = factory
        
        with fits.open(cube_file) as hdu_list:
            n_readouts = len(hdu_list[1].data['XBIN'][0])
            
            new_hdu_list = factory.construct_fits(template_file, n_samples=n_readouts, 
                                                  n_sp_kernels=3) # TODO: Need to dynamically set this.
            self.new_hdu_list = new_hdu_list
            
            #Populate Primary Header
            self.populate_primary_header(new_hdu_list[0], output_file, pds_label_file)
            
            channel = new_hdu_list[0].header['CHANNEL']
            
            # Get the Pandas data frame that contains the template information
            template = factory.get_template()
            template = template[~template['HDU NAME'].str.contains('PRIMARY')]
            hdu_names = template['HDU NAME']
            field_names = template['FIELD NAME']
            cg_field_names = template['CG FIELD NAME']
            
            # Copy data from the cube file to the new fits file.
            for cg_field_name, hdu_name, field_name in zip(cg_field_names, hdu_names, field_names):
                if type(cg_field_name) is not str:
                    continue
                
                new_data = new_hdu_list[hdu_name].data[field_name]
                
                # Insert cube data into the new fits file.  Dimensions will need to be reordered.
                if hdu_name == 'DATA':
                    for hdu_number in (1, 2, 3):
                        cube_data = np.moveaxis(hdu_list[hdu_number].data[cg_field_name], 
                                                (0, 1, 2, 3), (1, 0, 2, 3))
                        new_data[:, hdu_number-1, :, :] = np.squeeze(cube_data)
                elif (hdu_name == 'TIME') or (hdu_name == 'SC_GEOM'):
                    for hdu_number in (1, 2, 3):
                        if cg_field_name == 'UTC':
                            array_shape = new_data[:, hdu_number-1].shape
                            value = np.chararray(shape=array_shape, itemsize=26)
                            time = hdu_list[hdu_number].data[cg_field_name].squeeze()
                            for i in range(time.shape[0]):
                                value[i] = datetime.strptime(time[i], "%Y %b %d %H:%M:%S.%f").isoformat()
                        else:
                            value = np.moveaxis(hdu_list[hdu_number].data[cg_field_name],
                                                 (0,1), (1,0))
                        new_data[:, hdu_number-1] = value.squeeze()
                elif hdu_name == 'CAL':
                    cube_data = hdu_list[hdu_number].data[cg_field_name][0, 0, :1024, :].squeeze()
                    new_data[0, :, :] = cube_data
                    print(new_data)
                elif hdu_name == 'TARGET_GEOM':
                    for hdu_number in (1, 2, 3):
                        cube_shape = hdu_list[hdu_number].data[cg_field_name].shape
                        cube_data = hdu_list[hdu_number].data[cg_field_name].squeeze()
                        if len(cube_shape) == 2:
                            new_data[:, hdu_number-1] = cube_data
                        else:
                            new_data[:, hdu_number-1, :] = cube_data
                elif hdu_name == 'FOV_GEOM':
                    field_shape = new_data.shape
                    for hdu_number in (1, 2, 3):
                        if len(field_shape) == 4:
                            cg_corner_field_names = cg_field_name.split('|')
                            corner_index = 0
                            for cg_fn in cg_corner_field_names:
                                cube_data = hdu_list[hdu_number].data[cg_fn]
                                new_data[:, hdu_number-1, :, corner_index] = np.squeeze(cube_data)
                                corner_index += 1
                        else:
                            cube_data = hdu_list[hdu_number].data[cg_field_name].squeeze()
                            new_data[:, hdu_number-1, :] = np.squeeze(cube_data)
                else:
                    
                    # All values in the array should be the same.  Do a quick
                    # check of that here and raise an exception if it's not.
                    data_min = hdu_list[1].data[cg_field_name].min()
                    data_max = hdu_list[1].data[cg_field_name].max()
                    if data_min != data_max:
                        raise ValueError("All values in " + cg_field_name + "are not the same!")
                    
                    new_data[0] = data_max
                    
            # Add wavelengths, since they are not in the original cube
            xbin = hdu_list[1].data['XBIN'][0][0]
            table = Table(new_hdu_list['WAVELENGTH'].data)
            name = new_hdu_list['WAVELENGTH'].name
            header = new_hdu_list['WAVELENGTH'].header
            if channel == 'FUV':
                wavelengths = cassini_uvis_fuv_wavelengths(xbin)
                new_hdu_list['WAVELENGTH'].data['WAVELENGTH_FUV'] = wavelengths
                table = Table(new_hdu_list['WAVELENGTH'].data)
                table.remove_column('WAVELENGTH_EUV')
            else:
                wavelengths = cassini_uvis_euv_wavelengths(xbin)
                new_hdu_list['WAVELENGTH'].data['WAVELENGTH_EUV'] = wavelengths
                table = Table(new_hdu_list['WAVELENGTH'].data)
                table.remove_column('WAVELENGTH_FUV')
            new_hdu_list['WAVELENGTH'] = fits.BinTableHDU(table, name=name, header=header)
            
            # Kernel files
            if kernel_list is not None:
                self.populate_kernels_hdu(kernel_list)
                
        # Add obs_tick and obs_seconds information.
        self.set_obs_time(pds_label_file)
        
        # Set number of samples (readouts)
        self.set_number_of_samples()
        
    def set_obs_time(self, pds_label_file):
        '''
        Note that the sclk time is formatted in the following way:
        1/1609504792.123
        Where the "1" indicates the partition.  The number after the slash indicates
        the number of seconds.  And then after the period, it's the number of
        ticks, i.e. 1/256 of a second.
        '''
        pds_label = pvl.load(pds_label_file)
        sclk_time = pds_label['SPACECRAFT_CLOCK_START_COUNT'].split('/')[1]
        parts = sclk_time.split('.')
        seconds = parts[0]
        sub_seconds = parts[1]
        self.new_hdu_list['CONFIG'].data['OBS_SECONDS'] = seconds
        self.new_hdu_list['CONFIG'].data['OBS_TICKS'] = sub_seconds
    
    def set_number_of_samples(self):
        # Get the time from the header
        ephemeris_time = self.new_hdu_list['TIME'].data['TIME_ET']
        NZ = ephemeris_time.shape[0]
        self.new_hdu_list['CONFIG'].data['NUMBER_OF_SAMPLES'] = NZ
            
    def populate_kernels_hdu(self, kernel_list):
        '''
        Populate the spice kernels HDU.
        
        :param kernel_list:    List of spice kernels used for this observation
        :type kernel_list:     List of strings.
        '''
        kernels = np.array(kernel_list)
        
        map = {}
        map['LS_KRN'] = '.tls'
        map['SCL_KRN'] = '.tsc'
        map['SP_KRN'] = '.bsp'
        map['F_KRN'] = '.tf'
        map['PC_KRN'] = '.tpc'
        map['C_KRN'] = '.bc'
        map['INST_KRN'] = '.ti'
        
        table = Table(self.new_hdu_list['KERNELS'].data)
        name = self.new_hdu_list['KERNELS'].name
        header = self.new_hdu_list['KERNELS'].header
        
        # Scoop out everything but the comments, otherwise the old
        # dimensions will persist, creating bad file pointers.
        comments = [c for c in header.cards if c[0] == 'COMMENT']
        
        # Pop the placeholder Kernels table, so we can make a new one.
        self.new_hdu_list.pop('KERNELS')
        
        cols = []
        for key, value in map.items():
            
            w = np.flatnonzero(np.core.defchararray.find(kernels, value) != -1)
            
            # Use these kernels
            ks = kernels[w]
            
            # Get the max string length for this column
            max_len = max([len(kernel) for kernel in ks])
            
            # Set format
            data_format = 'A' + str(max_len)
            
            # Make the column objects
            col = fits.Column(name=key, format=data_format, array=ks[:], ascii=True)
            cols.append(col)
            
        # self.new_hdu_list['KERNELS'] = fits.BinTableHDU(table, name=name)
        self.new_hdu_list.insert(7, fits.TableHDU.from_columns(cols, name=name))
        
        # Add the comments back in
        for c in comments:
            self.new_hdu_list['KERNELS'].header['COMMENT'] = c[1]
    
    def write_to_file(self, output_file, overwrite=True):
        '''
        Write to output file.
        
        :param output_file:      Path to output FITS file.
        :type output_file:       pathLib.Path object.
        :param overwrite:        Set to overwrite output file.
        :type overwrite:         boolean
        '''
        self.factory.write_to_file(output_file, overwrite=overwrite)
       
    def populate_primary_header(self, hdu, cube_file, pds_label_file):
        '''
        Populate the primary header object.
        
        :param hdu:
        :type hdu:
        :param cube_file:
        :type cube_file:
        :param pds_label_file:
        :type pds_label_file:
        '''
        
        pds_label = pvl.load(pds_label_file)
        
        names = []
        values = []
        
        names.append('FILENAME')
        values.append(cube_file.name)
        
        names.append('PROD_ID')
        values.append(pds_label['PRODUCT_ID'])
        
        names.append('DATE')
        date_string = datetime.utcnow().isoformat()
        values.append(date_string)
        
        names.append('MISSION')
        values.append('Cassini')
        
        names.append('INSTRUME')
        values.append(pds_label['INSTRUMENT_NAME'])
        
        names.append('OBS_ID')
        values.append(pds_label['OBSERVATION_ID'])
        
        names.append('MPHASE')
        values.append(pds_label['MISSION_PHASE_NAME'])
        
        names.append('TRGTNAME')
        values.append(pds_label['TARGET_NAME'])
        
        names.append('OBS_UTC')
        utc = pds_label['START_TIME'].isoformat()
        values.append(utc)
        
        names.append('OBS_ET')
        cspice.furnsh('../kernels/lsk/naif0012.tls') # TODO: Dynamically download?
        t = cspice.str2et(utc.split('+')[0])
        values.append(t)
        cspice.unload('../kernels/lsk/naif0012.tls')
        
        names.append('END_UTC')
        values.append(pds_label['STOP_TIME'].isoformat())
        
        names.append('CHANNEL')
        values.append(pds_label['PRODUCT_ID'][:3])
        
        #TODO: Read version number from template or something.
        names.append('VERSION')
        values.append(1.0)
        
        sclk_time = np.double(pds_label['SPACECRAFT_CLOCK_START_COUNT'].split('/')[-1])
        orbit_number = get_orbit_number(sclk_time)
        names.append('ORBNUM')
        values.append(orbit_number)
        
        for name, value in zip(names, values):
            hdu.header[name] = value
            
            
class TitanCubeToFITS(CGCubeToFITS):
    '''
    This class maps an existing Cassini UVIS Cube file into the new FITS format
    for the Titan data products.
    '''
    
    def convert_to_fits(self, template_file, cube_file, output_file, pds_label_file, 
                        spice_dir,
                        detector_image_files=None, 
                        kernel_list=None):
        '''
        Convert a UVIS Cube file into the new FITS format.
        
        :param template_file:    Path to spreadsheet containing the template.
        :type template_file:     pathLib.Path object.
        :param cube_file:        Path to original cube file.
        :type cube_file:         pathLib.Path object.
        :param output_file:      Path to output FITS file.
        :type output_file:       pathLib.Path object.
        :param pds_label_file:   Path to PDS Label file for this observation
        :type pds_label_file:    pathLib.Path object.
        :param kernel_list:      List of kernel files used (for metadata purposes).
        :type kernel_list:       List of strings.
        
        '''
        
        # Call the superclass method to set up the FITS file general HDUs
        super().convert_to_fits(template_file, cube_file, output_file, pds_label_file, 
                                kernel_list=kernel_list)
        
        # Get HDUList defined in superclass.
        new_hdu_list = self.new_hdu_list
        
        # Based on the channel, exclude irrelevant fields.
        channel = new_hdu_list[0].header['CHANNEL']
        if channel == 'FUV':
            new_hdu_list.pop('DETECTOR_IMG_EUV')
        else:
            new_hdu_list.pop('DETECTOR_IMG_FUV')
        
        # Detector Images
        ymin = new_hdu_list['CONFIG'].data['IMG_YMIN'][0]
        ymax = new_hdu_list['CONFIG'].data['IMG_YMAX'][0]
        if detector_image_files is not None:
            for file in detector_image_files:
                # Read the data from the detector image
                detector_image = genfromtxt(file, delimiter=',')
                
                # Get the field name
                field_name = file.name[:-4].upper()
                
                # Write to the HDU
                hdu_name = 'DETECTOR_IMG_' + channel
                new_hdu_list[hdu_name].data[field_name][:, :ymax-ymin] = detector_image
                
        # Compute the SATURN_LOCAL_TIME value
        self.compute_saturn_local_time(spice_dir)
    
    def compute_saturn_local_time(self, spice_dir):
        
        # Get the time from the header
        ephemeris_time = self.new_hdu_list['TIME'].data['TIME_ET']
        
        # Get the kernels
        leap_second_kernel = spice_dir / 'lsk' / self.new_hdu_list['KERNELS'].data['LS_KRN'][0]
        planetary_constants_kernels = self.new_hdu_list['KERNELS'].data['PC_KRN']
        condition1 = np.core.defchararray.startswith(planetary_constants_kernels, 'pck')
        condition2 = np.core.defchararray.startswith(planetary_constants_kernels, 'cpck_rock')
        w = np.where(~(condition1 | condition2))
        planetary_constants_kernel = spice_dir / 'pck' / planetary_constants_kernels[w][0]
        ephemeris_kernels = self.new_hdu_list['KERNELS'].data['SP_KRN']
        condition = (np.core.defchararray.find(ephemeris_kernels, 'R_SCPSE_') != -1)
        w = np.where(condition)
        ephemeris_kernel = spice_dir / 'spk' / ephemeris_kernels[w][0]
        
        # Load the kernels
        cspice.furnsh(str(leap_second_kernel))
        cspice.furnsh(str(planetary_constants_kernel))
        cspice.furnsh(str(ephemeris_kernel))
        
        # Loop over all elements of ephemeris_time, convert to saturn local time.
        NT = ephemeris_time.shape[1]
        NZ = ephemeris_time.shape[0]
        for nt in range(NT):
            for nz in range(NZ):
                et = ephemeris_time[nz, nt]
                
                # psat_sun_iau is the vector from Saturn center to the Sun in the IAU_SATURN reference frame
                psat_sun_iau, sat_sun_ltime = cspice.spkpos('SUN', et, 
                                                            'IAU_SATURN', 'NONE', 'SATURN')
                
                # psat_ttn_iau is the vector from Saturn center to Titan in the IAU_SATURN reference frame
                psat_ttn_iau, sat_sun_ltime = cspice.spkpos('TITAN', et, 
                                                            'IAU_SATURN', 'NONE', 'SATURN')
                
                # Compute angles
                sun_angle = np.rad2deg(np.arctan2(psat_sun_iau[1], psat_sun_iau[0]))
                titan_angle = np.rad2deg(np.arctan2(psat_ttn_iau[1], psat_ttn_iau[0]))
                titan_sun_angle = ((titan_angle - sun_angle) + 180.00) % 360
                
                #Convert degree to hour
                slt_t = titan_sun_angle * 24.0 / 360.0
                
                # Copy the slt to the HDU
                self.new_hdu_list['TARGET_GEOM'].data['SATURN_LOCAL_TIME'][nz, nt] = slt_t
    
        # Unload kernels
        cspice.kclear()
    
def get_kernels_for_day(year, doy, spice_dir):
    '''
    Function to retrieve the spice kernels for a particular day.
    
    TODO: something less hacky.  The IDL-Python bridge is still broken on Mac
    as of IDL 8.8.1 so we'll have to just use subprocess.run to get the kernels.
    Better yet, just port the kernel_finder to Python.
    '''
    
    # Run IDL to get the list of kernels, dumping output to a temp directory
    idl_path_cmd = "!path = !path + ':" + str(Path.home() / "git" / "cassini-uvis-tools" / "kernel_finder") + "'"
    idl_kernel_cmd = "cassini_spice_kernel_list, " + str(year) + ", " + str(doy) + \
        ", fkernel_cg, fkernel, kernel_list, LOCAL_KERNEL_PATH='" + str(spice_dir) + "/'"
    idl_cmd = '/Applications/harris/envi56/idl88/bin/idl -e "' + idl_path_cmd + \
        ' & ' + idl_kernel_cmd + '"'
    cp = subprocess.run(idl_cmd, capture_output=True, shell=True)
    
    # Parse out the list of kernel files from the stdout.
    kernel_list = np.array(str(cp.stdout).split('\\n'))
    w = np.where('*****************************************************' == kernel_list)[0]
    kernel_list = kernel_list[w[0]+6:w[1]-1]
    
    
    # Call PlanetaryPy kernel finder.
    # year_doy = str(year) + '-' + str(doy)
    # kernels = list_kernels_for_day('cassini', year_doy, year_doy)
    # kernels = [kernel for kernel in kernels if not (kernel.endswith('.ti') and not 'uvis' in kernel)]
    
    
    # return the list of kernel files.  We don't need the files themselves,
    # only the list for the metadata.
    return kernel_list
    
if __name__ == '__main__':
    
    # orb_num = get_orbit_number(1871614519.013)
    # print(orb_num)
    
    spice_dir = Path('..') / 'spice'
    
    # kernel_list = get_kernels_for_day(2014, 265, spice_dir)
    # print(kernel_list)
    
    kernel_list = get_kernels_for_day(2005, 46, spice_dir)
    print(kernel_list)
    
    # Test list to avoid downloading for now.
    # kernel_list = ['naif0012.tls', 'cas00172.tsc', 'de432s.bsp',
    #      '150122R_SCPSE_14251_14283.bsp', 'cas_rocks_v18.tf', 'cas_status_v04.tf',
    #      'cas_v43.tf', '14212_14279py_as_flown.bc', '14263_14268ra.bc',
    #      '14244_14273ca_ISS.bc', 'cpck15Dec2017.tpc', 'cpck_rock_21Jan2011.tpc',
    #      'pck00010.tpc', 'cas_uvis_v07.ti']
    
    data_dir = Path('..') / 'data'
    template_dir = Path('..') / 'templates'
    template_file = template_dir / 'Titan_UVIS_data_definition_v0.2.xlsx'
    
    # cube_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined.fits'
    # new_fits_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined_new.fits'
    # label_file = data_dir / 'FUV2014_265_11_15.LBL'
    # detector_image_dir = data_dir / 'detector_images' / 'fuv'
    
    cube_file = data_dir / 'EUV2005_046_08_55_19_UVIS_003TI_EUVFUV001_PRIME_combined_CGFITS.fits'
    new_fits_file = data_dir / 'EUV2005_046_08_55_19_UVIS_003TI_EUVFUV001_PRIME_combined.fits'
    label_file = data_dir / 'EUV2005_046_08_55.LBL'
    detector_image_dir = data_dir / 'detector_images' / 'euv'
    
    detector_image_files = [f for f in detector_image_dir.iterdir()]
    
    
    
    converter = TitanCubeToFITS()
    converter.convert_to_fits(template_file, cube_file, new_fits_file, label_file, 
                              spice_dir,
                              detector_image_files=detector_image_files, 
                              kernel_list=kernel_list)
    converter.write_to_file(new_fits_file, overwrite=True)
    
    # Now, read the data to test
    # import matplotlib.pyplot as plt
    # import matplotlib.colors as colors
    # with fits.open(new_fits_file) as hdul:
    #
    #     xbin = hdul['CONFIG'].data['IMG_XBIN']
    #
    #     ymin = hdul['CONFIG'].data['IMG_YMIN']
    #     ymax = hdul['CONFIG'].data['IMG_YMAX']
    #
    #     # utc = hdul['TIME'].data['TIME_UTC']
    #     # print(utc)
    #     # et = hdul['TIME'].data['TIME_ET']
    #     # print(et)
    #
    #     for column in hdul['KERNELS'].data.columns:
    #         print(hdul['KERNELS'].data[column.name])
    #
    #     print(hdul['TARGET_GEOM'].data['SATURN_LOCAL_TIME'].shape)
    #     print(hdul['TARGET_GEOM'].data['SATURN_LOCAL_TIME'])
    #
    #     print("ymin = " + str(ymin))
    #     print("ymax = " + str(ymax))
    #
    #     # Get a the UVIS data
    #     data = np.squeeze(hdul['DATA'].data['UVIS_CALIBRATED'][:,0,:,:])
    #
    #     # Plot the image, summed in the NZ (samples) dimension
    #     dsum = np.sum(data, 0, where=(data > 0))
    #     plt.imshow(dsum, norm=colors.LogNorm(vmin=0.1, vmax=dsum.max()))
    #
    #     # Plot one of the detector images
    #     plt.figure()
    #     d = hdul['DETECTOR_IMG_FUV'].data['LYMAN_ALPHA']
    #     plt.imshow(d)
    #
    #     # Plot the image, summed along the wavelength dimension.
    #     plt.figure()
    #     dsum = np.sum(data, 1, where=(data > 0))
    #     plt.imshow(dsum, norm=colors.LogNorm(vmin=0.1, vmax=dsum.max()))
    #
    #     # Plot a spectrum, by summing over the NZ (sample) and NY (spatial) dimensions
    #     plt.figure()
    #     wavelengths = cassini_uvis_fuv_wavelengths(xbin)
    #     dsum = np.sum(data, 0, where=(data > 0))
    #     dsum = np.sum(dsum, 1, where=(dsum > 0))
    #     plt.plot(wavelengths, dsum)
    #
    #     # Plot a geometry array
    #     plt.figure()
    #     geo_data = hdul['FOV_GEOM'].data['RAYHEIGHT'][:,0,:].squeeze()
    #     plt.imshow(geo_data)
    #
    #     plt.show()
    #

