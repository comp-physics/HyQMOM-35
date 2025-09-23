clear
clc

syms M000 positive
syms umean vmean wmean real
syms C200 C110 C101 C020 C011 C002 real
syms C300 C210 C201 C120 C111 C102 C030 C021 C012 C003 real
syms C400 C310 C301 C220 C211 C202 C130 C121 C112 C103 C040 C031 C022 C013 C004 real

syms C4 [5 5 5] real

M4 = 0*C4;

%
C4(1,1,1) = 1;
%
C4(2,1,1) = 0;
C4(1,2,1) = 0;
C4(1,1,2) = 0;
%    
C4(3,1,1)=C200;
C4(2,2,1)=C110;
C4(2,1,2)=C101;
C4(1,3,1)=C020;
C4(1,2,2)=C011;
C4(1,1,3)=C002;
%
C4(4,1,1)=C300;
C4(3,2,1)=C210;
C4(3,1,2)=C201;
C4(2,3,1)=C120;
C4(2,2,2)=C111;
C4(2,1,3)=C102;
C4(1,4,1)=C030;
C4(1,3,2)=C021;
C4(1,2,3)=C012;
C4(1,1,4)=C003;
%
C4(5,1,1)=C400;
C4(4,2,1)=C310;
C4(4,1,2)=C301;
C4(3,3,1)=C220;
C4(3,2,2)=C211;
C4(3,1,3)=C202;
C4(2,4,1)=C130;
C4(2,3,2)=C121;
C4(2,2,3)=C112;
C4(2,1,4)=C103;
C4(1,5,1)=C040;
C4(1,4,2)=C031;
C4(1,3,3)=C022;
C4(1,2,4)=C013;
C4(1,1,5)=C004;

% find M4 from C4
for i=0:4 
    for j=0:4
        for k=0:4
            if i+j+k <= 4
                M4(i+1,j+1,k+1) = 0;
                for m = 0:i
                    for n = 0:j
                        for p = 0:k
                            M4(i+1,j+1,k+1) = M4(i+1,j+1,k+1) + M000*nchoosek(i,m)*nchoosek(j,n)*nchoosek(k,p)*(umean^(i-m))*(vmean^(j-n))*(wmean^(k-p))*C4(1+m,1+n,1+p);
                        end
                    end
                end
                simplify(M4(i+1,j+1,k+1));
            end
        end
    end
end

var = [M000,umean,vmean,wmean,C200,C110,C101,C020,C011,C002,...
               C300,C210,C201,C120,C111,C102,C030,C021,C012,C003,...
               C400,C310,C301,C220,C211,C202,C130,C121,C112,C103,C040,C031,C022,C013,C004];

matlabFunction(M4,'File','C4toM4_3D','Vars',var)