% To export current figure to image file
% Args (optional):
%   image_dump          name for image image dump file
%   image_dump_ext      file extension for image dump, in MATLAB format
%                       (https://uk.mathworks.com/help/matlab/ref/print.html)
function export_figure(varargin)
    image_dump_default = "plot";
    image_dump_ext_default = "-djpeg";

    ip = inputParser;
    addOptional(ip,'image_dump',image_dump_default,@isstring);
    addOptional(ip,'image_dump_ext',image_dump_ext_default,@isstring);

    parse(ip,varargin{:});
    file = ip.Results.image_dump;
    ext = ip.Results.image_dump_ext;

    fprintf("Dumping figure to file %s with MATLAB extension specifier %s", file, ext);
    print(file,ext);
end