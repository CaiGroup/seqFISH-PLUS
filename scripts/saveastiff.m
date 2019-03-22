function saveastiff(data, path, options)
% options.color
%   : true or FALSE
%   : If this is true, third dimension should be 3 and the data is saved as a color image.
% options.comp
%   : 'no', 'lzw', 'jpeg' or 'adobe'.
%     Compression type.
%       'no'    : Uncompressed(Default)
%       'lzw'   : lossless LZW
%       'jpeg'  : lossy JPEG
%       'adobe' : lossless Adobe-style
% options.message
%   : true or FALSE.
%     If this is false, all messages are skipped. 
% options.ask
%   : true or FALSE. If this is false, always overwrite an existing file.
% options.append
%   : true or FALSE
%     If path is exist, the data is appended to an existing file.
%     If path is not exist, this options is ignored.
% options.big 
%   : true or FALSE, 
%     Use 64 bit addressing and allows for files > 4GB
% 
% Defalut value of 'options' is
%     options.color   = false;
%     options.comp    = 'no';
%     options.message = false;
%     options.ask     = false;
%     options.append  = false;
%     options.big     = false;
% 
% Copyright (c) 2012, YoonOh Tak
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Gwangju Institute of Science and Technology (GIST), Republic of Korea nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

errcode = -1;
try
if isreal(data) == false
    errcode = 8; assert(false);
end
if nargin < 3
    options.color = false;
    options.comp = 'no';
    options.message = false;
    options.ask = false;
    options.append = false;
end
if isfield(options, 'message') == 0
    options.message = false;
end
if isfield(options, 'ask') == 0
    options.ask = false;
end
if isfield(options, 'append') == 0
    options.append = false;
end
if isfield(options, 'comp') == 0
    options.comp = 'no';
end
if isfield(options, 'color') == 0
    options.color = false;
end
if (options.color == false && ndims(data) > 3) ...
    || (options.color == true && ndims(data) > 4)
    errcode = 10; assert(false);
end
if isfield(options, 'big') == 0 
    options.big = false; 
end

[pathstr, fname, fext] = fileparts(which(path));
if ~isempty(pathstr)
    if strcmpi(pathstr, pwd)
        if options.append
            tappend = Tiff(path, 'r+');
        else
            if options.ask
                reply = input('File already exist, do you wish to overwrite? Y/n: ', 's');
                if isempty(reply), reply = 'Y'; end
                if any(upper(reply) ~= 'Y')
                    errcode = 1; assert(false);
                end
            end
        end
    else
        if options.append
            if options.ask
                display(['File path: ' pathstr '\' fname fext]);
                reply = input('File already exist but not in the current folder, do you wish to append? Y/n: ', 's');
                if isempty(reply), reply = 'Y'; end
                if any(upper(reply) == 'Y')
                    errcode = 1; assert(false);
                end
                tappend = Tiff(path, 'r+');
            end
        else
            if options.ask
                display(['File path: ' pathstr '\' fname fext]);
                reply = input('File already exist but not in the current folder, do you wish to overwrite? Y/n: ', 's');
                if isempty(reply), reply = 'Y'; end
                if any(upper(reply) ~= 'Y')
                    errcode = 1; assert(false);
                end
            end
        end
    end
else
    if exist(path, 'dir') && isempty(fname)
        errcode = 0; assert(false);
    end
    options.append = false;
end

if isempty(data)
    errcode = 2; assert(false);
end

if ~options.color
    if ndims(data) >= 4, errcode = 3; assert(false); end;
    [height, width, depth] = size(data);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
else
    if ndims(data) >= 5, errcode = 3; assert(false); end;
    [height, width, rgb, depth] = size(data);
    if rgb ~= 3, errcode = 4; assert(false); end;
    if isa(data, 'uint8'), data = uint8(data); end;
    tagstruct.Photometric = Tiff.Photometric.RGB;
end

tagstruct.ImageLength = height;
tagstruct.ImageWidth = width;
tagstruct.SamplesPerPixel = (options.color+1)*(options.color+2)/2;
tagstruct.RowsPerStrip = height;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;

% Compresstion type : http://en.wikipedia.org/wiki/Tagged_Image_File_Format
switch lower(options.comp)
    case 'no'
        tagstruct.Compression = Tiff.Compression.None;
    case 'lzw'
        tagstruct.Compression = Tiff.Compression.LZW;
    case 'jpeg'
        tagstruct.Compression = Tiff.Compression.JPEG;
    case 'adobe'
        tagstruct.Compression = Tiff.Compression.AdobeDeflate;
    otherwise
        if options.ask
            reply = input('Unsupported compression type : Data will not be compressed.\nDo you wish to continue? Y/n: ', 's');
            if isempty(reply), reply = 'Y'; end
            if any(upper(reply) ~= 'Y')
                tagstruct.Compression = Tiff.Compression.None;
            else
                errcode = 0; assert(false);
            end
        else
            tagstruct.Compression = Tiff.Compression.None;
        end
end

switch class(data)
    case {'uint8', 'uint16', 'uint32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'int16', 'int32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        if options.color == 3
            errcode = 8; assert(false);
        end
    case {'single', 'double'}
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        errcode = 9; assert(false);
end

switch class(data)
    case {'uint8', 'int8'}
        tagstruct.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tagstruct.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tagstruct.BitsPerSample = 32;
    case {'single'}
        tagstruct.BitsPerSample = 32;
    case {'double'}
        tagstruct.BitsPerSample = 64;
    otherwise
        errcode = 9; assert(false);
end

tStart = tic;
if ~options.append
    if ispc
        tempfile = '~$temporal.tif';
    else
        tempfile = '._temporal.tif';
    end
    s=whos('data'); 
    if s.bytes > 2^32-1 || options.big 
        t = Tiff(tempfile, 'w8'); 
    else 
        t = Tiff(tempfile, 'w'); 
    end
    if ispc
        fileattrib(tempfile, '+h +w', '', 's');
    end
    for d = 1:depth
        t.setTag(tagstruct);
        if ~options.color
            t.write(data(:, :, d));
        else
            t.write(data(:, :, :, d));
        end
        if d ~= depth
            t.writeDirectory();
        end
    end

    t.close();
    if ispc
        fileattrib(tempfile, '-h +w');
    end
    if isempty(pathstr)
        movefile(tempfile, path, 'f');
    else
        movefile(tempfile, [pathstr '\' fname fext], 'f');
    end
else
    while ~tappend.lastDirectory();
        tappend.nextDirectory();
    end
    
    try
        tappend.writeDirectory();
        for d = 1:depth
            tappend.setTag(tagstruct);
            if ~options.color
                tappend.write(data(:, :, d));
            else
                tappend.write(data(:, :, :, d));
            end
            if d ~= depth
                tappend.writeDirectory();
            end
        end
    catch
        errcode = 11; assert(false);
    end
    tappend.close();
end

tElapsed = toc(tStart);
if options.message
    display(sprintf('File saved successfully. Elapsed time : %.3f s.', tElapsed));
end

catch exception
    if exist('t', 'var')
        t.close();
        if exist('tappend', 'var'), tappend.close(); end
        delete(tempfile);
    end
    
    if options.message
        switch errcode
            case 0
                error 'Invalide path.';
            case 1
            case 2
                error '''data'' is empty.';
            case 3
                error 'Data dimension is too large.';
            case 4
                error 'Third dimesion (color depth) should be 3.';
            case 6
                error 'RGB color image can not have int8, int16 or int32 format.';
            case 8
                error 'It does not support complex numbers.';
            case 9
                error 'Unsupported data type.';
            case 10
                error 'Dimension of source data is too large.';
            case 11
                error 'If file size exceeded, create a new big tiff file using ''options.big = true'''
            otherwise
                rethrow(exception);
        end
    end
end
end
