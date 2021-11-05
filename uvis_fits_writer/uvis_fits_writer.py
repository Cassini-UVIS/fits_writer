'''
Created on 1 Nov 2021

@author: Josh Elliott : joshua.elliott@lasp.colorado.edu
'''

from astropy.io import fits
from astropy.io import ascii
from astropy.time import Time
from astropy.table import Table
import numpy as np

class UVISFitsHDUList(fits.HDUList):
    '''
    This class defines an HDUList object.
    '''
    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        
        # Create an empty primary HDU
        self.primary_hdu = fits.PrimaryHDU()
        self.append(self.primary_hdu)
        
    @property
    def get_primary_hdu(self):
        return self.primary_hdu
        
    def add_header_item(self, hdu, name, value, comment=None):
        header = hdu.header
        header[name] = (value, comment)
        
    def __add_new_data(self,  hdu_class, **kwargs):
        hdu = hdu_class(**kwargs)
        self.append(hdu)
        return hdu
    
    def add_new_image(self, **kwargs):
        return self.__add_new_data(fits.ImageHDU, **kwargs)
        
    def add_new_binary_table(self, **kwargs):
        return self.__add_new_data(fits.BinTableHDU, **kwargs)

class UVISFitsTemplateFactory(object):
    '''
    This class defines an object which writes a FITS formatted data file
    given a data structure definition spreadsheet.
    '''
    
    def __init__(self):
        pass
        
    def construct_fits_file(self, fits_file_name, csv_template_file):
        # Parse the template
        csv_template = ascii.read(csv_template_file, encoding='utf-8-sig')
        
        # Instantiate an HDUList, which contains an empty PrimaryHDU
        hdu_list = UVISFitsHDUList()
        
        # Construct the primary header
        self.construct_primary_header(hdu_list, csv_template)
        
        # Create additional HDUs based on the template CSV file
        self.construct_hdus(hdu_list, csv_template)
        
        # Write to file
        hdu_list.writeto(fits_file_name, overwrite=True)
    
    def construct_primary_header(self, hdu_list, csv_template):
        '''
        Method to construct header information to an HDU
        
        Parameters
        ----------
        hdu : PrimaryHDU
            HDU to write header information to
        csv_template : astropy.table
            Template row definitions.
        Returns
        -------
        '''
        
        primary_hdu = hdu_list.get_primary_hdu
        
        w = np.where(csv_template['HDU NAME'] == 'PRIMARY_HEADER')
        header_rows = csv_template[w]
        for row in header_rows:
            dtype = row['DATA TYPE']
            if dtype == 'bytes':
                dtype = np.str
                value = ''
            else:
                value = np.array([0], dtype=dtype)[0]
            comment = row['COMMENT'] # TODO: Some comments are too long, need to break up?
            name = row['FIELD NAME']
            hdu_list.add_header_item(primary_hdu, name, value, comment=comment)
    
    def construct_hdu(self, hdu_list, hdu_name, hdu_rows):
        w = np.where(hdu_rows['HDU NAME'] == hdu_name)
        rows = hdu_rows[w]
        
        data = {}
        for row in rows:
            field_name = row['FIELD NAME']
            dtype = row['DATA TYPE']
                
            # Get the dimensions
            dimensions = row['DIMENSION']
            column_length = row['COLUMN LENGTH']
            dimensions = (dimensions + ',' + column_length).split(',')
            
            array_shape = []
            for dimension in dimensions:
                # TODO: Make this more generic.  For now use standard UVIS dimensions.
                if dimension == 'NX':
                    array_shape.append(1024)
                elif dimension == 'NY':
                    array_shape.append(64)
                elif dimension == 'NZ':
                    array_shape.append(1) # For the template, just do a single readout.
                elif dimension == 'NT':
                    array_shape.append(3)
                else:
                    # Numeric type
                    array_shape.append(int(dimension))
            
            # Data types.
            if dtype == 'bytes':
                dtype = np.str #TODO: use dtype=object here?
                
            value = np.zeros(shape=array_shape, dtype=dtype)
            
            if value.size == 1:
                value = value[0]
                
            data[field_name] = value
            
            # TODO: Comments into the header.
            #comment = row['COMMENT'] 
        
        try:
            data = Table(data)
        except:
            print('bad')
        hdu_list.add_new_binary_table(data=data, name=hdu_name)
    
    def construct_hdus(self, hdu_list, csv_template):
        # Get HDU names
        w = np.where((csv_template['HDU NAME'] != 'PRIMARY_HEADER') & 
                     (csv_template['HDU NAME'] != 'PRIMARY'))
        
        hdu_rows = csv_template[w]
        hdu_names = np.unique(hdu_rows['HDU NAME'])
        
        for name in hdu_names:
            self.construct_hdu(hdu_list, name, hdu_rows)
        
    
if __name__ == '__main__':
    csv_file = 'Titan_UVIS_data_definition_v0.1.csv'
    new_fits_file = 'test_fits_file.fits'
    factory = UVISFitsTemplateFactory()
    factory.construct_fits_file(new_fits_file, csv_file)
    
    hdul = fits.open(new_fits_file)
    hdul.info()
    