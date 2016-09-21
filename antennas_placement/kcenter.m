function [S, R, XY] = kcenter(exname, limit)
%KCENTER - algorithm to place the two types of antennas
%
% Syntax:  [S, R, XY] = kcenter(exname, limit)
%
% Inputs:
%    exname - excel file
%    limit  - 1300 meters for MICRO antennas
%           -  200 meters for PICO antennas
%
% Outputs:
%   S - indices of antennas
%   R - covering radius
%   XY possitions of antennas
%
% Example: 
%   [S R XY] = kcenter('buildingsIDXY.xlsx', 1300)
%   [S R XY] = kcenter('buildingsPointsPart1.xlsx', 200)
%   [S R XY] = kcenter('buildingsPointsPart2.xlsx', 200)
% 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: none

% Author        : Marios Kyriakou
% Work address  : KIOS Research Center, University of Cyprus
% email         : mkiria01@ucy.ac.cy
% Website       : http://www.kios.ucy.ac.cy
% Last revision : September 2016

%------------- BEGIN CODE --------------

[~,~,res] = xlsread(exname);
[~,ex]=fileparts(exname);
% save results
file=[ex,'Res.csv']; 

ID = res(:,1);
X1 = str2num(char(res{:,2}));
Y1 = str2num(char(res{:,3}));
P = [X1 Y1];
S=1;

% greedy optimisation
n = length(P);
L = zeros(1,n);
L(1) = S;
% initial DD, R
Skns = zeros(n, length(P));
R = zeros(n-1,1);

% calculate symmetric matrix A
A=dist(P,P');
Skns(1, :) = A(1,:);   
SknsMin = Skns(1, :);

for i=2:n %868
  [R(i-1), newL] = max(SknsMin);
  L(i) = newL;
  Skns(i,:) = A(L(i),:);
  SknsMin = min(SknsMin, Skns(i,:));
end

%results
R(n) = max(SknsMin);
S =L';
XY = P(L,:);
S=S(find(R>limit));

% write CSV file
if exist(file)
   delete(file) 
end

csvwrite(file,[S,XY(find(R>limit),1),XY(find(R>limit),2)]);
[tlines]=regexp( fileread(file), '\n', 'split');
tlinesNew = cell(1,length(tlines)+1);
tlinesNew{1}='ID, X, Y';
tlinesNew(2:end)=tlines;
f = fopen(file,'w');
for i=1:length(tlinesNew)
    fprintf(f,tlinesNew{i});
    fprintf(f,'\n');
end
fclose all;

%------------- END OF CODE --------------
