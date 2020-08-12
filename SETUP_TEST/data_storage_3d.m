function path2data = data_storage_3d()
% Usage: path2data = data_storage_3d()
% 
% Purpose: Provides path to data storage location on computer that
%          currently runs the code (edit this file to add new computers).
%
% Input:
%   none
%
% Output:
%   path2data : [char] : path to data storage
%
% Part of M3TET - 3D FINITE ELEMENT MANTLE CONVECTION CODE
%
% Developed by J.Hasenclever & J.Phipps Morgan, 2007-2010
% Email contact: joerg.hasenclever@zmaw.de, jhasenclever@geomar.de
% For numerical methods see online Ph.D. thesis of J.Hasenclever
% (http://www.sub.uni-hamburg.de/opus/volltexte/2010/4873/)
%
% JH Jan 2012
% JH Feb 2016 : added new machines
%

hostname = get_hostname;
switch hostname
    case 'b3pc14'
        path2data = 'D:/jhasenclever/';
    case 'b3pc16'
        path2data = 'G:/jhasenclever/';
    case 'b3pc17'
        path2data = 'D:/Projects/';
    case 'Maupiti2'
        path2data = 'J:/Projects/';
    case 'fuego'
        path2data = '/Data_Fu/jhasenclever/';
    case 'colima'
        path2data = '/Data_Fu/jhasenclever/';
    case 'Snaefell'
        path2data = '/Users/joha/TMP/';
    case 'jorgeMAC'
        path2data = '/Users/jorge/Tests/';
    case 'jorgeLAPTOP'
        path2data = 'C:\Users\Tara\Documents\PhD\Tests\';
    case 'clusterRHUL'
        path2data = '/mlbraid/jorge';
    case 'jmc1'
        path2data = '/mlbraid/jorge';
    case 'clusterRHUL_2'
        path2data = '/mldata2/jorge';
    case 'glmsc12'
        path2data = '/mldata2/jorge';
    otherwise
        if strcmp(hostname(1:4),'node') || strcmp(hostname,'master')
            path2data = './OUTPUT/';
        elseif strcmp(hostname(1:4),'rzcl')
            path2data = './OUTPUT/';
        elseif strcmp(hostname(1:4),'nesh')
            path2data = './OUTPUT/';
        else
            error(' No location for output data defined for this computer.\n You have to edit file "data_storage_3d.m"!');
        end
end

% '\' can cause problems on unix system; however matlab can handle '/' on
% all architectures
path2data(path2data=='\') = '/';
% make sure there's a slash at the end of the path
if path2data(end)~='/'
    path2data = [path2data '/'];
end

end % END OF FUNCTION data_storage_3d