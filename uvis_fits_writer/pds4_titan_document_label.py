'''
Created on Sep 20, 2022

@author: Josh Elliott : joshua.elliott@lasp.colorado.edu
'''

from astropy.io import fits
from xml.etree import ElementTree as ET
from uvis_fits_writer.pds4_label_base import PDS4Label
from pathlib import Path
from datetime import datetime

def remove_namespace(xmldoc, namespace):
    ns = u'{%s}' % namespace
    nsl = len(ns)
    for elem in xmldoc.getiterator():
        if elem.tag.startswith(ns):
            elem.tag = elem.tag[nsl:]

class PDS4TitanDocumentLabelCreator(PDS4Label):
    '''
    This class defines a PDS4 label for a document.
    '''
    
    def __init__(self, xml_template, fits_file, document_file, output_dir=Path('.')):
        '''
        Constructor
        '''
        
        # Call superclass constructor.
        super().__init__()
        
        # Get the hdu list and primary hdu references.
        self.fits_file = fits_file
        self.hdu_list = fits.open(self.fits_file)
        self.primary_hdu = self.hdu_list['PRIMARY']
        
        self.document_file = document_file
        
        if 'plot' in xml_template.name:
            self.doc_type = 'plot'
            if 'eps' in xml_template.name:
                self.file_extention = '.eps'
            if 'png' in xml_template.name:
                self.file_extention = '.png'
        else: # animation
            self.doc_type = 'animation'
            self.file_extention = '.mp4'
        
        # Get the preamble ready
        self.get_xml_preamble()
        
        # Set the namespace
        ns = 'http://pds.nasa.gov/pds4/pds/v1'
        ET.register_namespace('', ns)
                              
        self.xml_root = ET.parse(xml_template).getroot()
        
    def create_pds4_label(self):
        
        # Loop over each element in the XML tree
        for elem in self.xml_root.iter():
            # Modify the template with the relevant document information.
            if elem.tag.endswith('logical_identifier'):
                # Here we need to append the file name to the ID.
                parts = elem.text.split(':')
                parts[-1] = self.document_file.name[:-4].lower().replace('-', '_')
                elem.text = ':'.join(parts)
            if elem.tag.endswith('title'):
                # Modify the title.
                elem.text = elem.text + self.fits_file.name[:-5].replace('-', '_')
            if elem.tag.endswith('Citation_Information'):
                # Citation information.  Use the current year in UTC.
                for elem2 in elem.iter():
                    if elem2.tag.endswith('publication_year'):
                        elem2.text = str(datetime.utcnow().year)
                    if elem2.tag.endswith('description'):
                        elem2.text = elem2.text + self.fits_file.name[:-5].replace('-', '_')
            if elem.tag.endswith('Modification_History'):
                # Mod history.  In this case, use the current year-month-day.
                for elem2 in elem.iter():
                    if elem2.tag.endswith('modification_date'):
                        elem2.text = \
                            str(datetime.utcnow().year) + "-" + \
                            str(datetime.utcnow().month) + "-" + \
                            str(datetime.utcnow().day) + "Z"
            if elem.tag.endswith('Context_Area'):
                for elem2 in elem.iter():
                    if elem2.tag.endswith('start_date_time'):
                        elem2.text = self.format_time_stamp(self.primary_hdu.header['OBS_UTC'])
                    if elem2.tag.endswith('stop_date_time'):
                        elem2.text = self.format_time_stamp(self.primary_hdu.header['END_UTC'])
            if elem.tag.endswith('Reference_List'):
                for elem2 in elem.iter():
                    if elem2.tag.endswith('lid_reference'):
                        elem2.text = elem2.text + self.fits_file.name[:-5].lower().replace('-', '_')
            if elem.tag.endswith('}Document'): # The } is necessary here so it doesn't confuse with 'Product Document'
                for elem2 in elem.iter():
                    if elem2.tag.endswith('document_name'):
                        elem2.text = elem2.text + self.fits_file.name[:-5].replace('-', '_')
                    if elem2.tag.endswith('copyright'):
                        elem2.text = str(datetime.utcnow().year)
                    if elem2.tag.endswith('publication_date'):
                        elem2.text = elem2.text = \
                            str(datetime.utcnow().year) + "-" + \
                            str(datetime.utcnow().month) + "-" + \
                            str(datetime.utcnow().day) + "Z"
                    if elem2.tag.endswith('description'):
                        elem2.text = elem2.text + self.fits_file.name[:-5].replace('-', '_')
                    if elem2.tag.endswith('file_name'):
                        elem2.text = self.document_file.name[:-4] + self.file_extention
        
        # root = self.xml_root.getroot()
        # elt = root.find('//{' + ns + '}Product_Document')
        # print(elt)
        # ET.dump(elt)
        
        # root = self.xml_root.getroot()
        # xmlstr = ET.tostring(root, encoding='utf8', method='xml')
        # print(xmlstr)
        ET.dump(self.xml_root)
        
if __name__ == '__main__':
    data_dir = Path('..') / 'data'
    output_dir = Path('..') / 'output'
    if not output_dir.exists():
        output_dir.mkdir()
    template_dir = Path('..') / 'templates'
    xml_template = template_dir / 'titan_plot_template_eps.xml'
    fits_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined_new.fits'
    label_file = output_dir / ('documentxml_' + fits_file.name[:-5] + '_plot.xml')
    
    # eps
    document_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined_plot.eps'
    label_creator = PDS4TitanDocumentLabelCreator(xml_template, fits_file, document_file)
    label_creator.create_pds4_label()
    label_creator.write_to_file(document_file.name[:-4] + '.xml')
    
    # png
    xml_template = template_dir / 'titan_plot_template_png.xml'
    label_file = output_dir / ('documentxml_' + fits_file.name[:-5] + '.xml')
    document_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined_animation_frame_000.png'
    label_creator = PDS4TitanDocumentLabelCreator(xml_template, fits_file, document_file)
    label_creator.create_pds4_label()
    label_creator.write_to_file(document_file.name[:-4] + '.xml')
    
    # mp4
    document_file = data_dir / 'FUV2014_265_11_15_21_UVIS_208TI_EUVFUV002_PRIME_combined_animation.mp4'
    xml_template = template_dir / 'titan_animation_template.xml'
    label_file = output_dir / ('documentxml_' + fits_file.name[:-5] + '_animation.xml')
    label_creator = PDS4TitanDocumentLabelCreator(xml_template, fits_file, document_file)
    label_creator.create_pds4_label()
    label_creator.write_to_file(document_file.name[:-4] + '.xml')
        