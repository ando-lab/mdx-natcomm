classdef PilatusHeader
    % ref http://www.bernstein-plus-sons.com/software/CBF/doc/cif_img_1.7.10.html
    % and https://www.dectris.com/technical_pilatus.html?file=tl_files/root/support/technical_notes/pilatus/Pilatus_CBF_Header_Specification.pdf
    properties(SetAccess = protected)
        filepath
        date_time
        HeaderContents = struct([]);
        CIFHeader = struct([]);
        binary_start
    end
    properties(Access = protected)
        cif_lines = {}
        header_lines = {}
        non_binary_length = 4096
        rawheader
    end
    properties(Constant, Access = protected)
        SPACE_EQUIVALENT_CHARACTERS = get_space_equivalent_characters();
        CIF_BINARY_KEYWORDS = get_cif_binary_keywords();
        NON_OPTIONAL_KEYWORDS = get_non_optional_keywords();
        OPTIONAL_KEYWORDS = get_optional_keywords();
    end
    methods (Static = true) % a useful function for reading header contents
        function h = read(varargin)
            hobj = io.PilatusHeader(varargin{:});
            h = hobj.HeaderContents;
        end
    end
    methods
        function obj = PilatusHeader(filepath,non_binary_length)
            assert(nargin>=1,'at least one argument expected');
            obj.filepath = filepath;
            if nargin==2
                obj.non_binary_length = non_binary_length;
            end
            
            obj.rawheader = obj.read_header_data();
            obj.binary_start = obj.seek_binary_start();
            obj.cif_lines = obj.read_cif_lines();
            obj.CIFHeader = obj.parse_cif_header();
            
            if obj.has_pilatus_cbf_convention()
                obj.header_lines = obj.read_header_lines();
                obj.HeaderContents = obj.parse_header();
                obj.date_time = obj.get_date_time();
            end
        end
        
        value = read_header_data(obj)

        tf = has_pilatus_cbf_convention(obj)
        
        value = seek_binary_start(obj)
        
        value = read_cif_lines(obj)
        
        value = read_header_lines(obj)

        value = spaced_header_lines(obj)

        Header = parse_header(obj)
        
        value = datatype_handling(obj,values,keyword,datatype)

        value = get_date_time(obj)
        
        Header = parse_cif_header(obj)
 
    end
end