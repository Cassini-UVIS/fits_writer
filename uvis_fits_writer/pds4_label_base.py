'''
Created on Sep 20, 2022

@author: Josh Elliott : joshua.elliott@lasp.colorado.edu
'''

from pathlib import Path
import xml.etree.ElementTree as ET
from xml.dom import minidom
from abc import ABC, abstractmethod # Abstract Base Class

class PDS4Label(ABC):
    '''
    This class defines a base class for a PDS4 label.
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
    
    def __init__(self, product_type='Product_Observational'):
        '''
        Constructor
        '''
        
        self.schema_version = '1I00'
        self.schema_disp_ver = '1I00_1510'
        self.information_model_version = '1.16.0.0'
        self.product_type = product_type
        
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
        self.xml_root = ET.Element(self.product_type, 
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
    
    @abstractmethod
    def create_pds4_label(self):
        raise NotImplementedError
    
    def format_time_stamp(self, time_stamp):
        '''
        The PDS requires the time stamp to be in ISO8601 format, UTC time zone.
        TODO: check on the number of decimal places allowed.
        '''
        return time_stamp[:24] + 'Z'
    
    def write_to_file(self, output_file):      
        # https://stackoverflow.com/a/14493981/2313806
        def pretty_print(data):
            return '\r\n'.join([line for line in minidom.parseString(data).toprettyxml(indent='  ').split('\n') if line.strip()])
          
        # Prettify the XML string
        xmlstr = ET.tostring(self.xml_root, encoding='utf-8', xml_declaration=None)
        #xmlstr = minidom.parseString(xmlstr).toprettyxml(indent='\t')
        xmlstr = pretty_print(xmlstr)
        
        # Replace generic preamble with the one we want.
        # TODO: There has got to be a better way to do this:
        xmlstr = xmlstr.replace('<?xml version="1.0" ?>', self.preable_content)
        
        with open(output_file, 'w') as f:
            f.write(xmlstr)
    