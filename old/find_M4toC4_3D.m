clear
clc

syms M000 positive
syms M100 M010 M001 real
syms M200 M110 M101 M020 M011 M002 real
syms M300 M210 M201 M120 M111 M102 M030 M021 M012 M003 real
syms M400 M310 M301 M220 M211 M202 M130 M121 M112 M103 M040 M031 M022 M013 M004 real

syms M4 [5 5 5] real

M4 = 0*M4;
C4 = 0*M4;

umean=M100/M000;
vmean=M010/M000;
wmean=M001/M000;

%
M4(1,1,1)=M000;
%
M4(2,1,1)=M100;
M4(1,2,1)=M010;
M4(1,1,2)=M001;
%    
M4(3,1,1)=M200;
M4(2,2,1)=M110;
M4(2,1,2)=M101;
M4(1,3,1)=M020;
M4(1,2,2)=M011;
M4(1,1,3)=M002;
%
M4(4,1,1)=M300;
M4(3,2,1)=M210;
M4(3,1,2)=M201;
M4(2,3,1)=M120;
M4(2,2,2)=M111;
M4(2,1,3)=M102;
M4(1,4,1)=M030;
M4(1,3,2)=M021;
M4(1,2,3)=M012;
M4(1,1,4)=M003;
%
M4(5,1,1)=M400;
M4(4,2,1)=M310;
M4(4,1,2)=M301;
M4(3,3,1)=M220;
M4(3,2,2)=M211;
M4(3,1,3)=M202;
M4(2,4,1)=M130;
M4(2,3,2)=M121;
M4(2,2,3)=M112;
M4(2,1,4)=M103;
M4(1,5,1)=M040;
M4(1,4,2)=M031;
M4(1,3,3)=M022;
M4(1,2,4)=M013;
M4(1,1,5)=M004;

% find C4 from M4
for i=0:4 
    for j=0:4
        for k=0:4
            if i+j+k <= 4
                C4(i+1,j+1,k+1) = 0;
                for m = 0:i
                    for n = 0:j
                        for p = 0:k
                            C4(i+1,j+1,k+1) = C4(i+1,j+1,k+1) + nchoosek(i,m)*nchoosek(j,n)*nchoosek(k,p)*((-umean)^(i-m))*((-vmean)^(j-n))*((-wmean)^(k-p))*M4(1+m,1+n,1+p)/M000;
                        end
                    end
                end
                simplify(C4(i+1,j+1,k+1));
            end
        end
    end
end

var = [M000,M100,M200,M300,M400,M010,M110,M210,M310,M020,M120,M220,M030,M130,M040,...
       M001,M101,M201,M301,M002,M102,M202,M003,M103,M004,M011,M111,M211,M021,M121,...
       M031,M012,M112,M013,M022];

matlabFunction(C4,'File','M4toC4_3D','Vars',var)