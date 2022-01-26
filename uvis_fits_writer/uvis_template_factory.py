'''
Created on 1 Nov 2021

@author: Josh Elliott : joshua.elliott@lasp.colorado.edu

'''

from astropy.io import fits
from astropy.io import ascii
from astropy.table import Table
import numpy as np
from _collections import OrderedDict
import pandas
from pathlib import Path

class UVISFitsHDUList(fits.HDUList):
    '''
    This class defines an HDUList object.
    '''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        # Create an empty primary HDU
        self.primary_hdu = fits.PrimaryHDU()
        self.append(self.primary_hdu)
        
    def get_primary_hdu(self):
        return self.primary_hdu
        
    def add_header_item(self, hdu, name, value, comment=None):
        hdu.header[name] = (value, comment)
        
    def __add_new_data(self,  hdu_class, **kwargs):
        hdu = hdu_class(**kwargs)
        self.append(hdu)
        return hdu
    
    def add_new_image(self, **kwargs):
        return self.__add_new_data(fits.ImageHDU, **kwargs)
        
    def add_new_binary_table(self, **kwargs):
        return self.__add_new_data(fits.BinTableHDU, **kwargs)
    
    def add_new_ascii_table(self, **kwargs):
        return self.__add_new_data(fits.TableHDU, **kwargs)

class UVISFitsTemplateFactory(object):
    '''
    This class defines an object which writes a FITS formatted data file
    given a data structure definition spreadsheet.
    '''
    
    def __init__(self):
        pass
        
    def construct_fits(self, template_file, n_samples=1, n_sp_kernels=1):
        # Parse the template
        self.template = pandas.read_excel(str(template_file))
        self.template.replace(np.nan, '') # Prefer empty strings to nans
        
        # Instantiate an HDUList, which contains an empty PrimaryHDU
        self.hdu_list = UVISFitsHDUList()
        
        # Construct the primary header
        self.construct_primary_header(self.template)
        
        # Create additional HDUs based on the template file
        self.construct_hdus(self.template, n_samples=n_samples, n_sp_kernels=n_sp_kernels)
        
        return self.hdu_list
    
    def get_template(self):
        return self.template
        
    def write_to_file(self, fits_file_name, overwrite=True):
        # Write to file
        self.hdu_list.writeto(fits_file_name, overwrite=overwrite)
    
    def construct_primary_header(self, template):
        '''
        Method to construct header information to an HDU
        
        Parameters
        ----------
        template : astropy.table
            Template row definitions.
        Returns
        -------
        '''
        
        primary_hdu = self.hdu_list.get_primary_hdu()
        
        # w = np.where(template['HDU NAME'] == 'PRIMARY_HEADER')
        # header_rows = template[w]
        header_rows = template[template['HDU NAME'] == 'PRIMARY_HEADER']
        for row in header_rows.index:
            dtype = header_rows['DATA TYPE'][row]
            if dtype == 'bytes':
                dtype = np.str
                value = ''
            else:
                value = np.array([0], dtype=dtype)[0]
            comment = header_rows['COMMENT'][row] # TODO: Some comments are too long, need to break up?
            name = header_rows['FIELD NAME'][row]
            self.hdu_list.add_header_item(primary_hdu, name, value, comment=comment)
    
    def construct_hdu(self, hdu_name, hdu_rows, n_samples=1, n_sp_kernels=1):
        rows = hdu_rows[hdu_rows['HDU NAME'] == hdu_name]
        
        header = fits.Header()
        
        # Create a dict to store the data
        data = OrderedDict()
        
        # Loop over rows
        for row in rows.index:
            field_name = rows['FIELD NAME'][row]
            dtype = rows['DATA TYPE'][row]
                
            # Get the dimensions
            dimensions = rows['DIMENSION'][row]
            column_length = rows['COLUMN LENGTH'][row]
            dimensions = (str(column_length) + ',' + str(dimensions)).split(',')
            
            array_shape = []
            for dimension in dimensions:
                # TODO: Make this more generic.  For now use standard UVIS dimensions.
                # TODO: How to handle dimensions of 1?  Should we np.squeze() them?
                if dimension == 'NX':
                    array_shape.append(1024)
                elif dimension == 'NY':
                    array_shape.append(64)
                elif dimension == 'NZ':
                    array_shape.append(n_samples)
                elif dimension == 'NT':
                    array_shape.append(3)
                elif dimension == 'NK':
                    array_shape.append(n_sp_kernels)
                else:
                    array_shape.append(int(dimension))
                    
            if hdu_name == 'WAVELENGTH':
                array_shape = array_shape[1:]
                    
            if array_shape[-1] == 1:
                array_shape = array_shape[0:-1]
            
            if field_name == 'TIME_UTC':
                value = np.chararray(shape=array_shape, itemsize=26)
            elif hdu_name == 'KERNELS':
                value = np.chararray(shape=(1, n_sp_kernels), itemsize=80)
            else:
                # Data types.
                if dtype == 'bytes':
                    dtype = np.str #TODO: use dtype=object here?
                value = np.zeros(shape=array_shape, dtype=dtype)
                
            data[field_name] = value
            
            # TODO: Find a way to add this to the TTYPE* lines instead
            header['COMMENT'] = field_name + ' : ' + rows['COMMENT'][row]
  
        # Convert to binary table object.
        data = Table(data, masked=True)
        
        # TODO: For Kernels, perhaps use ASCII Table instead?
        self.hdu_list.add_new_binary_table(data=data, name=hdu_name, header=header)
    
    def construct_hdus(self, template, n_samples=1, n_sp_kernels=1):
        # Get HDU names
        hdu_rows = template[
            (template['HDU NAME'] != 'PRIMARY_HEADER') & 
            (template['HDU NAME'] != 'PRIMARY')
            ]
        hdu_names = np.unique(hdu_rows['HDU NAME'])
        
        for name in hdu_names:
            self.construct_hdu(name, hdu_rows, n_samples=n_samples, n_sp_kernels=n_sp_kernels)
    
if __name__ == '__main__':
    template_dir = Path('..') / 'templates'
    template_file = template_dir / 'Titan_UVIS_data_definition_v0.2.xlsx'
    new_fits_file = 'test_fits_file.fits'
    factory = UVISFitsTemplateFactory()
    factory.construct_fits(template_file, n_samples=20)
    factory.write_to_file(new_fits_file, overwrite=True)
    
    hdul = fits.open(new_fits_file)
    hdul.info()
    