'''
Created on 10 Mar 2022

@author: Josh Elliott : joshua.elliott@lasp.colorado.edu

'''

from astropy.io import fits
from pathlib import Path
import numpy as np
import xml.etree.ElementTree as ET
from xml.dom import minidom
import csv
import pandas

class PDS4LabelCreator(object):
    '''
    This class defines a PDS4 label given a FITS file data product. 
    '''
    
    # Note: The FITS standard specifies MSB, i.e. big-endian, for ALL data types.
    type_map = {
        # Unsigned int types
        'uint8':   {'name':'UnsignedByte', 'size':1, 'offset':0},
        'uint16':  {'name':'UnsignedMSB2', 'size':2, 'offset':0},
        'uint32':  {'name':'UnsignedMSB4', 'size':4, 'offset':0},
        'uint64':  {'name':'UnsignedMSB8', 'size':8, 'offset':0},
        
        # Signed int types 
        'int8':    {'name':'SignedByte', 'size':1, 'offset':0},
        'int16':   {'name':'SignedMSB2', 'size':2, 'offset':int("8000", 16)},
        'int32':   {'name':'SignedMSB4', 'size':4, 'offset':int("80000000", 16)},
        'int64':   {'name':'SignedMSB8', 'size':8, 'offset':int("8000000000000000", 16)},
        
        # Float types
        'float32':  {'name':'IEEE754MSBSingle', 'size':4, 'offset':0},
        'float64': {'name':'IEEE754MSBDouble', 'size':8, 'offset':0},
        
        # String types
        'string': {'name':'ASCII_String', 'size':1, 'offset':0}
        }
    
    def __init__(self, fits_file, template_file, data_product_level='1A'):
        '''
        Constructor
        '''
        
        # TODO: Check which schema version numbers we should be using.
        # These versions listed below are the most current on the PDS website:
        # https://pds.nasa.gov/datastandards/dictionaries/index-1.18.0.0.shtml
        self.schema_version = '1I00'
        self.schema_disp_ver = '1I00_1510'
        self.information_model_version = '1.16.0.0'
        
        # Data Product Level
        # TODO: Check with team on the proper data product levels.
        self.product_level_number = data_product_level
        self.pds_product_level = self.get_pds_product_level(data_product_level)
        
        # Get the hdu list and primary hdu references.
        self.fits_file = fits_file
        self.hdu_list = fits.open(self.fits_file)
        self.primary_hdu = self.hdu_list['PRIMARY']
        
        # Get the preamble ready
        self.get_xml_preamble()
        
        # Get the template information:
        self.template = pandas.read_excel(str(template_file), 
                                          sheet_name=['L2A_product_definition', 'HDU Descriptions'])
        print('test')
    
    def get_xml_preamble(self):        
        # Read the XML preamble:
        # TODO: There has to be a better way to add these via the ElementTree library
        #     At present all I can find is a way to add the first line via the
        #     xml_declaration=True keyword.  But the others I'm not sure.
        preamble_file = Path('..') / 'pds' / 'xml_preamble.xml'
        with open(preamble_file, 'r') as f:
            self.preable_content = f.read()
    
    def create_xml_root(self): # pds4_xml_model:
        # TODO: Collections, Bundles, etc.
        self.xml_root = ET.Element('Product_Observational', 
                {
                    'xmlns': 'http://pds.nasa.gov/pds4/pds/v1',
                    'xmlns:pds': 'http://pds.nasa.gov/pds4/pds/v1',
                    'xmlns:disp': 'http://pds.nasa.gov/pds4/disp/v1',
                    'xmlns:xsi': 'http://www.w3.org/2001/XMLSchema-instance',
                    'xsi:schemaLocation': 'http://pds.nasa.gov/pds4/pds/v1      ' + \
                        'http://pds.nasa.gov/pds4/pds/v1/PDS4_PDS_1G00.xsd     ' + \
                        'https://pds.nasa.gov/pds4/disp/v1      ' + \
                        'https://pds.nasa.gov/pds4/disp/v1/PDS4_DISP_1G00_1500.xsd'
                }
            )
        
    def get_pds_product_level(self, level, pds=False):
        '''
        Data product level definitions:
        https://pdsmgmt.gsfc.nasa.gov/documents/PDS_Policy_on_Data_Processing_Levels_20130311.pdf
        https://www.earthdata.nasa.gov/engage/open-data-services-and-software/data-information-policy/data-levels
        '''
        if level == '1A':
            product_level = 'Raw'
        elif level == '1B':
            product_level = 'Partially_Processed'
        elif level == '2':
            product_level = 'Calibrated'
        elif level == '3':
            product_level = 'Derived'
            
        return product_level
    
    def get_modification_history(self, xml_parent):
        # Mod history group
        mod_history = ET.SubElement(xml_parent, 'Modification_History')
        
        # Point to mod history file.
        mod_history_file = Path('..') / 'pds' / 'modification_history.txt'
        
        # Loop over each row in the CSV file and append to the mod history.
        with open(mod_history_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                mod_detail = ET.SubElement(mod_history, 'Modification_Detail')
                mod_date = ET.SubElement(mod_detail, 'modification_date')
                mod_date.text = row['date']
                mod_version = ET.SubElement(mod_detail, 'version_id')
                mod_version.text = row['version']
                mod_description = ET.SubElement(mod_detail, 'description')
                mod_description.text = row['description']
    
    def get_instrument(self):
        '''
        Get the instrument name.  Usually the short name, i.e. UVIS.
        TODO: Need a better way to get this.  Perhaps we should include this in
            the fits header.
        '''
        instrument = self.primary_hdu.header.cards['INSTRUME'].value
        if self.primary_hdu.header.cards['INSTRUME'].value == 'ULTRAVIOLET IMAGING SPECTROGRAPH':
            instrument = 'UVIS'
        return instrument
    
    def get_spacecraft(self):
        return self.primary_hdu.header.cards['MISSION'].value
    
    def get_target_name(self):
        return self.primary_hdu.header.cards['TRGTNAME'].value
        
    def get_product_lid(self):
        '''
        Create a product lid as described in section 5 of the PDS Data Provider's Handbook
        TODO: Need a more generic way to construct these.
        '''
        lid = 'urn:nasa:pds:cassini-uvis_titan-library:data:' + self.fits_file.stem.lower()
        return lid
    
    def get_investigation_area_lid(self):
        '''
        Create a lid for the investigation area.
        TODO: Need a more generic way to construct these.
        '''
        lid = 'urn:nasa:pds:context:investigation:mission.cassini-huygens'
        return lid
    
    def get_observing_system_instrument_lid(self):
        '''
        Create a lid for the observing system instrument.
        TODO: Need a more generic way to construct these.
        '''
        lid = 'urn:nasa:pds:context:instrument:uvis.co'
        return lid
    
    def get_observing_system_spacecraft_lid(self):
        '''
        Create a lid for the observing system spacecraft.
        TODO: Need a more generic way to construct these.
        '''
        lid = 'urn:nasa:pds:context:instrument_host:spacecraft.co'
        return lid
        
    def get_target_lid(self):
        '''
        Create a lid for the observing system.
        TODO: Need a more generic way to construct these.
        '''
        lid = 'urn:nasa:pds:context:target:satellite.saturn.titan'
        return lid
        
    def format_time_stamp(self, time_stamp):
        '''
        The PDS requires the time stamp to be in ISO8601 format, UTC time zone.
        TODO: check on the number of decimal places allowed.
        '''
        return time_stamp[:24] + 'Z'
    
    def create_sub_element(self, parent, name, keys=None, text=None):
        if keys is not None:
            child = ET.SubElement(parent, name, keys)
        else:
            child = ET.SubElement(parent, name)
        if text is not None:
            child.text = text
        return child
    
    def create_pds4_label(self, output_file):

        # Get the size of the primary HDU image, if present.
        have_primary = self.primary_hdu.data is not None
        if have_primary:
            primary_size = self.primary_hdu.data.shape
        
        # Create the root level of the XML file
        self.create_xml_root()
        
        # Create <Identification_Area> tag -------------------------------------
        identification_area = self.create_sub_element(self.xml_root, 'Identification_Area')
        
        logical_identifier = self.create_sub_element(identification_area, 'logical_identifier', 
                                           text=self.get_product_lid())
        
        version_id = self.create_sub_element(identification_area, 'version_id', 
                                             text = str(self.primary_hdu.header.cards['VERSION'].value))
        
        title = self.create_sub_element(identification_area, 'title', 
                                        text = self.primary_hdu.header.cards['MISSION'].value + ' ' + \
                                        self.primary_hdu.header.cards['INSTRUME'].value)
                     
        information_model_version = self.create_sub_element(identification_area, 'information_model_version', 
                                                            text = self.information_model_version)
        
        product_class = self.create_sub_element(identification_area, 'product_class', 
                                                text = self.xml_root.tag)
        
        self.get_modification_history(identification_area)
        
        # End <Identification_Area> tag ----------------------------------------
        
        # Create <Observation_Area> tag ----------------------------------------
        observation_area = self.create_sub_element(self.xml_root, 'Observation_Area')
        
        time_coordinates = self.create_sub_element(observation_area, 'Time_Coordinates')
        start_date_time = self.create_sub_element(time_coordinates, 'start_date_time',
                                                  text = self.format_time_stamp(self.primary_hdu.header['OBS_UTC']))
        stop_date_time = self.create_sub_element(time_coordinates, 'stop_date_time', 
                                                 text=self.format_time_stamp(self.primary_hdu.header['END_UTC']))
        
        primary_result_summary = self.create_sub_element(observation_area, 'Primary_Result_Summary')
        purpose = self.create_sub_element(primary_result_summary, 'purpose', 
                                          text = 'Science')
        processing_level = self.create_sub_element(primary_result_summary, 'processing_level', 
                                                   text=self.pds_product_level)
        
        science_facets = self.create_sub_element(primary_result_summary,'Science_Facets')
        wavelength_range = self.create_sub_element(science_facets, 'wavelength_range', 
                                                   text='Ultraviolet')
        domain = self.create_sub_element(science_facets, 'domain', 
                                         text = 'Atmosphere') # TODO: Need to be able to set this for other types of observations.
        discipline_name = self.create_sub_element(science_facets, 'discipline_name', 
                                                  text='Atmospheres') # TODO: Need to be able to set this for other types of observations.
        
        facet1 = self.create_sub_element(science_facets, 'facet1', text='Structure')
        
        investigation_area = self.create_sub_element(observation_area, 'Investigation_Area')
        name = self.create_sub_element(investigation_area, 'name', text='IUVS')
        type_ = self.create_sub_element(investigation_area, 'type', text='Mission')
        
        internal_reference = self.create_sub_element(investigation_area, 'Internal_Reference')
        lid_reference = self.create_sub_element(internal_reference, 'lid_reference', 
                                                text=self.get_investigation_area_lid())
        reference_type = self.create_sub_element(internal_reference, 'reference_type', 
                                                 text='data_to_investigation')
        
        observing_system = self.create_sub_element(observation_area, 'Observing_System')
        name = self.create_sub_element(observing_system, 'name', 
                                       text=self.get_spacecraft())
        observing_system_component = self.create_sub_element(observing_system, 'Observing_System_Component') # Instrument
        name = self.create_sub_element(observing_system_component, 'name', 
                                       text=self.get_instrument())
        type = self.create_sub_element(observing_system_component, 'type', 
                                       text='Instrument')
        internal_reference = self.create_sub_element(observing_system_component, 'Internal_Reference')
        lid_reference = self.create_sub_element(internal_reference, 'lid_reference', 
                                                text=self.get_observing_system_instrument_lid())
        reference_type = self.create_sub_element(internal_reference, 'reference_type',
                                                 text='is_instrument')
        observing_system_component = self.create_sub_element(observing_system, 'Observing_System_Component') # Spacecraft
        name = self.create_sub_element(observing_system_component, 'name', 
                                       text=self.get_spacecraft())
        type = self.create_sub_element(observing_system_component, 'type', 
                                       text='Spacecraft')
        internal_reference = self.create_sub_element(observing_system_component, 'Internal_Reference')
        lid_reference = self.create_sub_element(internal_reference, 'lid_reference', 
                                                text=self.get_observing_system_spacecraft_lid())
        reference_type = self.create_sub_element(internal_reference, 'reference_type',
                                                 text='is_instrument_host')
        
        target_identification = self.create_sub_element(observation_area, 'Target_Identification')
        name = self.create_sub_element(target_identification, 'name', 
                                       text=self.get_target_name())
        type = self.create_sub_element(target_identification, 'type', 
                                       text='TODO') # TODO: We need a way to get the target type.
        internal_reference = self.create_sub_element(target_identification, 'Internal_Reference')
        lid_reference = self.create_sub_element(internal_reference, 'lid_reference', 
                                                text=self.get_target_lid())
        reference_type = self.create_sub_element(internal_reference, 'reference_type',
                                                 text='data_to_target')
        
        # TODO: May need to add information here
        discipline_area = self.create_sub_element(observation_area, 'Discipline_Area')
        
        # End <Observation_Area> tag -------------------------------------------
        
        # Begin <File_Area_Observational> tag ----------------------------------
        file_area_observational = self.create_sub_element(self.xml_root, 'File_Area_Observational')
        file_grp = self.create_sub_element(file_area_observational, 'File')
        self.create_sub_element(file_grp, 'file_name', 
                                text=self.fits_file.name)
        self.create_sub_element(file_grp, 'local_identifier', 
                                text='file')
        self.create_sub_element(file_grp, 'creation_date_time', 
                                text=self.format_time_stamp(self.primary_hdu.header.cards['DATE'].value))
        
        # Loop over HDUs
        i = 0
        for hdu in self.hdu_list:
            # TODO: Primary HDU case.
            if i < 1:
                i += 1
                continue # skip the primary hdu
            self.create_table_element(file_area_observational, hdu)
            
        # End <File_Area_Observational> tag ------------------------------------
      
    def create_table_element(self, parent, hdu):
        # Get the byte offset of the hdu header
        file_info = hdu.fileinfo()
        header_offset = file_info['hdrLoc']
        data_offset = file_info['datLoc']
        header_size = data_offset - header_offset
        
        data = hdu.data
        
        fields = hdu.columns.names
        n_rows = data[fields[0]].shape[0] # The first dimension is the number of rows
        n_fields = len(fields)
        
        # total number of fields - sum of number of elements in each field
        # total record size - sum of bytes of each field
        # number of groups - number of fields which are themselves arrays
        n_fields_lbl = 0
        record_length = 0
        n_groups = 0
        field_pos = np.zeros((n_fields,), dtype=np.uint64)
        field_size = np.zeros((n_fields,), dtype=np.uint64)
        for i in range(n_fields):
            d = data[fields[i]][0]
            if hdu.name == 'DATA':
                print("shape of " + fields[i], d.shape)
            this_n_fields = d.size
            if d.dtype.name.startswith('str'):
                size = np.max(np.char.str_len(data[fields[i]]))
                if size < 1:
                    size = 1
                field_size[i] = size
            else:
                field_size[i] = self.type_map[d.dtype.name]['size']
            
            field_size[i] *= this_n_fields
            field_pos[i] = record_length+1
            record_length += field_size[i]
            if len(d.shape) > 0:
                n_groups += 1
            else:
                n_fields_lbl += 1
                
        # Get descriptions
        hdu_names = self.template['HDU Descriptions']['HDU NAME'].values
        rows = self.template['HDU Descriptions'][hdu_names == hdu.name]
        header_description = rows['Header'].values[0]
        table_description = rows['Table'].values[0]
        
        # Make the header for this table
        header_grp = self.create_sub_element(parent, 'Header')
        self.create_sub_element(header_grp, 'name', 
                                text=hdu.name + '_HEADER')
        self.create_sub_element(header_grp, 'offset', 
                                keys={'unit':'byte'}, 
                                text=str(header_offset))
        self.create_sub_element(header_grp, 'object_length', 
                                keys={'unit':'byte'}, 
                                text=str(header_size))
        self.create_sub_element(header_grp, 'parsing_standard_id', 
                                text='FITS 3.0')
        self.create_sub_element(header_grp, 'description', 
                                text=header_description)
        
        # Table
        table_binary_grp = self.create_sub_element(parent, 'Table_Binary')
        self.create_sub_element(table_binary_grp, 'name', 
                                text=hdu.name + '_TABLE')
        self.create_sub_element(table_binary_grp, 'offset', 
                                keys={'unit': 'byte'}, 
                                text=str(data_offset))
        self.create_sub_element(table_binary_grp, 'records', 
                                text=str(n_rows))
        self.create_sub_element(table_binary_grp, 'description',
                                text=table_description)
        
        # Table Record
        record_binary_grp = self.create_sub_element(table_binary_grp, 'Record_Binary')
        self.create_sub_element(record_binary_grp, 'fields', 
                                text=str(n_fields_lbl))
        self.create_sub_element(record_binary_grp, 'groups', 
                                text=str(n_groups))
        self.create_sub_element(record_binary_grp, 'record_length', 
                                text='%d' % record_length, 
                                keys={'unit': 'byte'})
        
        # Loop over each field in the HDU and add <Field_Binary> tags.
        record_parent = record_binary_grp
        for i_field in range(n_fields):
            d = data[fields[i_field]][0]
            if d.size > 1:
                d = d.squeeze()
            
            s = d.shape
            if len(s) == 0:
                s = (1,)
            
            # Hard case for arrays. We have a (possibly nested set of, or possibly zero) <Group_Field_Binary> 
            # tags describing each dimension of the array. The outermost tag corresponds to the last index,
            # and the innermost tag to the first.
    
            if d.dtype.name.startswith('str'):
                size = np.max(np.char.str_len(d))
                if size < 1:
                    size = 1
                this_fieldsize = size
            else:
                this_fieldsize = self.type_map[d.dtype.name]['size']
            
            for i_dim_pds in range(len(s)):         # Outermost group has i_dim_pds=0
                i_dim_idl = len(s) - 1 - i_dim_pds  # Innermost group has i_dim_idl=0
                
                # print nested group headers
                this_fields = (i_dim_idl == 0) # only the innermost dimension has a field. All the others have a group
                this_groups = 0 if this_fields else np.prod(s[0:i_dim_idl+1]) # TODO: double check the +1
                this_group_length = this_fieldsize * np.prod(s[0:i_dim_idl+1])
        
                if (s[i_dim_idl] <= 1) and not d.dtype.name.startswith('str') : # If we are a scalar or singleton array dimension, don't print the group header
                    field_location = field_pos[i_field]
                else: #is scalar else not scalar
                    # Print a group descriptor
                    field_location = 1
                    group_field_binary = self.create_sub_element(record_parent, 'Group_Field_Binary')
                    repetitions = self.create_sub_element(group_field_binary, 'repetitions', text=str(s[i_dim_idl]))
                    group_location = field_pos[i_field] if i_dim_pds == 0 else 1
                    n_group_fields = 1 if i_dim_idl == 0 else 0
                    n_group_groups = 1 - n_group_fields
                    fields_elt = self.create_sub_element(group_field_binary, 'fields', text=str(n_group_fields))
                    groups_elt = self.create_sub_element(group_field_binary, 'groups', text=str(n_group_groups))
                    
                    group_location_elt = self.create_sub_element(group_field_binary, 'group_location', 
                                                             text=str(group_location), 
                                                             keys={'unit': 'byte'})
                    
                    group_length = np.prod(s[0:i_dim_idl+1]) * this_fieldsize
                    group_length_elt = self.create_sub_element(group_field_binary, 'group_length', 
                                                             text=str(group_length), 
                                                             keys={'unit': 'byte'})
                    
                    # Set the new parent
                    record_parent = group_field_binary
                    
                if i_dim_idl == 0: # we are at the center
                    if d.dtype.name.startswith('str'):
                        this_offset = self.type_map['string']['offset']
                        this_data_type = self.type_map['string']['name']
                    else:
                        this_offset = self.type_map[d.dtype.name]['offset']
                        this_data_type = self.type_map[d.dtype.name]['name']
                    
                    field_binary_elt = self.create_sub_element(record_parent, 'Field_Binary')
                    name_elt = self.create_sub_element(field_binary_elt, 'name', text=fields[i_field])
                    field_location_elt = self.create_sub_element(field_binary_elt, 'field_location', text=str(field_location), keys={'unit': 'byte'})
                    data_type_elt = self.create_sub_element(field_binary_elt, 'data_type', text=this_data_type)
                    field_length_elt = self.create_sub_element(field_binary_elt, 'field_length', text=str(this_fieldsize), keys={'unit': 'byte'})
                    
                    
                    if this_offset != 0:
                        value_offset_elt = self.create_sub_element(field_binary_elt, 'value_offset', text=str(this_offset))
                        
            # Set the new parent
            record_parent = record_binary_grp
    
    def write_to_file(self):        
        # Prettify the XML string
        xmlstr = ET.tostring(self.xml_root, encoding='utf-8', xml_declaration=None)
        xmlstr = minidom.parseString(xmlstr).toprettyxml(indent='\t')
        
        # Replace generic preamble with the one we want.
        # TODO: There has got to be a better way to do this:
        xmlstr = xmlstr.replace('<?xml version="1.0" ?>', self.preable_content)
        
        with open(label_file, 'w') as f:
            f.write(xmlstr)
    
if __name__ == '__main__':
    
    data_dir = Path('..') / 'data'
    fits_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined_new.fits'
    label_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined_new.xml'
    template_dir = Path('..') / 'templates'
    template_file = template_dir / 'Titan_UVIS_data_definition_v0.2.xlsx'
    
    label_creator = PDS4LabelCreator(fits_file, template_file)
    label_creator.create_pds4_label(label_file)
    label_creator.write_to_file()
    
    hdu_list = fits.open(fits_file)
    hdu_list.info()
    for hdu in hdu_list:
        print(hdu.name)
        print(hdu.filebytes())
        print(hdu.fileinfo())
    
    print(hdu_list['KERNELS'].data['INST_KRN'][0])
    
        