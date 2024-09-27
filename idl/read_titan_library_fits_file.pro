;+
; :Description:
;    Read a Cassini UVIS Titan library fits file from the following PDS archive:
;    https://atmos.nmsu.edu/data_and_services/atmospheres_data/TITAN/titan-lib.html
;
; :Params:
;    filename [in, required, String] : File path to the fits file.
;    
; :Keywords:
;    HEADER [out, optional, String[]] : String array containing the primary header information.
;    
; :Returns:
;   Data structure
;
; :Requires:
;   Written with IDL 9.0, but should run with anything > 8.5
;   IDLAstro : https://github.com/wlandsman/IDLAstro
;   
; :Author: 
;   Josh Elliott : joshua.p.elliott@jpl.nasa.gov
; 
; :History:
;   Created 27 Sept 2024
;-

function read_titan_library_fits_file, filename, HEADER=primary_header
  compile_opt idl2
  
  ; Get the number of hdus
  fits_info, filename, N_EXT=n_hdus, /SILENT

  ; Get the header
  !null = mrdfits(filename, 0, primary_header, /SILENT)

  hdus = hash()
  
  for i=1, n_hdus do begin ; start from 1, 0 is the "primary" which is just the header in our case.
    hdu = mrdfits(filename, i, header, /SILENT)
    w = where(header.startswith('EXTNAME'))
    hdu_name = strtrim((header[w].split("'"))[1], 2)
    hdu = orderedhash(hdu)
    hdus[hdu_name] = hdu.tostruct()
  endfor

  ; Return the struct
  return, hdus.tostruct()
end
