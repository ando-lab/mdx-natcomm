%% Gaussian Fits to Scattering Factors
%
% * Steve Meisburger
% * May 1, 2017
%
% Scattering factors (coherent and incoherent) were transcribed from the
% tables in:
%
% * Hubbell et al. Atomic form factors, incoherent scattering functions,
% and photon scattering cross sections, 1975. Journal of Physical and
% Chemical Reference Data 4, 471?538. doi:10.1063/1.555523
%
% They are fit using a four-gaussian function, as described in:
%
% * Brown PJ, Fox AG, Maslen EN, O'Keefe MA, Willis BTM, 2006. Intensity of
% diffracted intensities, in: Prince, E. (Ed.), International Tables for
% Crystallography Volume C: Mathematical, Physical and Chemical Tables.
% Springer, pp. 554?595. doi:10.1107/97809553602060000600
%
% On May 11, I noted that the coherent scattering factors for H did not
% match the values produced by the CCP4 program SFALL. This is because
% SFALL (correctly) uses the values for bonded hydrogen, which are
% different from free H2. The values can be found in table IV of Hubbel et
% al. (copied from earlier work by Setwart 1965 and Stewart and Bentley
% 1973). 
%% plot all available data in the tables
figure(1);clf
fn = dir('tables/*.dat');
for j=1:length(fn)
    d = dlmread(fullfile('tables',fn(j).name));
    plot(d(:,1),d(:,2),'.-','displayName',fn(j).name);hold on;
end
lg = legend('show');
set(lg,'Interpreter','none');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
%% Fit parameters
ngauss = 4; % fit to four gaussians plus constant
maxAngle = 2; % max value to fit, of sin(th)/lambda in A^{-1}
%%
% *Hydrogen*
Z = 1;
d = dlmread('tables/coh_H.dat');
coeffs_coh_H = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_H),'-');

d = dlmread('tables/incoh_H.dat');
coeffs_incoh_H = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_H(1:2:end) = -coeffs_incoh_H(1:2:end);
coeffs_incoh_H(end) = coeffs_incoh_H(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_H),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%%
% *Carbon*
Z = 6;
d = dlmread('tables/coh_C.dat');
coeffs_coh_C = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_C),'-');

d = dlmread('tables/incoh_C.dat');
coeffs_incoh_C = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_C(1:2:end) = -coeffs_incoh_C(1:2:end);
coeffs_incoh_C(end) = coeffs_incoh_C(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_C),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%%
% *Nitrogen*
Z = 7;
d = dlmread('tables/coh_N.dat');
coeffs_coh_N = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_N),'-');

d = dlmread('tables/incoh_N.dat');
coeffs_incoh_N = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_N(1:2:end) = -coeffs_incoh_N(1:2:end);
coeffs_incoh_N(end) = coeffs_incoh_N(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_N),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%%
% *Oxygen*
Z = 8;
d = dlmread('tables/coh_O.dat');
coeffs_coh_O = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_O),'-');

d = dlmread('tables/incoh_O.dat');
coeffs_incoh_O = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_O(1:2:end) = -coeffs_incoh_O(1:2:end);
coeffs_incoh_O(end) = coeffs_incoh_O(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_O),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%%
% *Phosphorus*
Z = 15;
d = dlmread('tables/coh_P.dat');
coeffs_coh_P = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_P),'-');

d = dlmread('tables/incoh_P.dat');
coeffs_incoh_P = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_P(1:2:end) = -coeffs_incoh_P(1:2:end);
coeffs_incoh_P(end) = coeffs_incoh_P(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_P),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%%
% *Sulfur*
Z = 16;
d = dlmread('tables/coh_S.dat');
coeffs_coh_S = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_S),'-');

d = dlmread('tables/incoh_S.dat');
coeffs_incoh_S = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_S(1:2:end) = -coeffs_incoh_S(1:2:end);
coeffs_incoh_S(end) = coeffs_incoh_S(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_S),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);

%%
% *Sodium*
Z = 11;
d = dlmread('tables/coh_Na.dat');
coeffs_coh_Na = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_Na),'-');

d = dlmread('tables/incoh_Na.dat');
coeffs_incoh_Na = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_Na(1:2:end) = -coeffs_incoh_Na(1:2:end);
coeffs_incoh_Na(end) = coeffs_incoh_Na(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_Na),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);

%%
% *Chlorine*
Z = 17;
d = dlmread('tables/coh_Cl.dat');
coeffs_coh_Cl = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_Cl),'-');

d = dlmread('tables/incoh_Cl.dat');
coeffs_incoh_Cl = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_Cl(1:2:end) = -coeffs_incoh_Cl(1:2:end);
coeffs_incoh_Cl(end) = coeffs_incoh_Cl(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_Cl),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);

%%
% *Argon*
Z = 18;
d = dlmread('tables/coh_Ar.dat');
coeffs_coh_Ar = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_Ar),'-');

d = dlmread('tables/incoh_Ar.dat');
coeffs_incoh_Ar = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_Ar(1:2:end) = -coeffs_incoh_Ar(1:2:end);
coeffs_incoh_Ar(end) = coeffs_incoh_Ar(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_Ar),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);

%%
% *Potassium*
Z = 19;
d = dlmread('tables/coh_K.dat');
coeffs_coh_K = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_K),'-');

d = dlmread('tables/incoh_K.dat');
coeffs_incoh_K = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_K(1:2:end) = -coeffs_incoh_K(1:2:end);
coeffs_incoh_K(end) = coeffs_incoh_K(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_K),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);

%%
% *Calcium*
Z = 20;
d = dlmread('tables/coh_Ca.dat');
coeffs_coh_Ca = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_Ca),'-');

d = dlmread('tables/incoh_Ca.dat');
coeffs_incoh_Ca = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_Ca(1:2:end) = -coeffs_incoh_Ca(1:2:end);
coeffs_incoh_Ca(end) = coeffs_incoh_Ca(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_Ca),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);

%%
% *Iron*
Z = 26;
d = dlmread('tables/coh_Fe.dat');
coeffs_coh_Fe = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_Fe),'-');

d = dlmread('tables/incoh_Fe.dat');
coeffs_incoh_Fe = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_Fe(1:2:end) = -coeffs_incoh_Fe(1:2:end);
coeffs_incoh_Fe(end) = coeffs_incoh_Fe(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_Fe),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%%
% *Cobalt*
Z = 27;
d = dlmread('tables/coh_Co.dat');
coeffs_coh_Co = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_Co),'-');

d = dlmread('tables/incoh_Co.dat');
coeffs_incoh_Co = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_Co(1:2:end) = -coeffs_incoh_Co(1:2:end);
coeffs_incoh_Co(end) = coeffs_incoh_Co(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_Co),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%%
% *Cadmium*
Z = 48;
d = dlmread('tables/coh_Cd.dat');
coeffs_coh_Cd = fitSumOfGaussians(d(:,1),d(:,2),ngauss,true,maxAngle);

figure(1);clf;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_coh_Cd),'-');

d = dlmread('tables/incoh_Cd.dat');
coeffs_incoh_Cd = fitSumOfGaussians(d(:,1),Z-d(:,2),ngauss,true,maxAngle);
coeffs_incoh_Cd(1:2:end) = -coeffs_incoh_Cd(1:2:end);
coeffs_incoh_Cd(end) = coeffs_incoh_Cd(end)+Z;

hold on;
plot(d(:,1),d(:,2),'o',d(:,1),sumOfGaussians(d(:,1),coeffs_incoh_Cd),'-');
xlabel('sin(\theta)/\lambda');ylabel('F(q,Z) or S(q,Z)');
title(['Z = ',num2str(Z)]);
%% print coefficients to table
fprintf(1,'\n%3s %5s ','','');
for j = 1:ngauss, fprintf(1,'%9s ',['a' num2str(j)],['b' num2str(j)]); end
fprintf(1,'%9s ','c');
fprintf(1,'\n%3s %5s ','H','coh');   fprintf(1,'%9.5f ',coeffs_coh_H);
fprintf(1,'\n%3s %5s ','H','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_H);
fprintf(1,'\n%3s %5s ','C','coh');   fprintf(1,'%9.5f ',coeffs_coh_C);
fprintf(1,'\n%3s %5s ','C','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_C);
fprintf(1,'\n%3s %5s ','N','coh');   fprintf(1,'%9.5f ',coeffs_coh_N);
fprintf(1,'\n%3s %5s ','N','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_N);
fprintf(1,'\n%3s %5s ','O','coh');   fprintf(1,'%9.5f ',coeffs_coh_O);
fprintf(1,'\n%3s %5s ','O','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_O);
fprintf(1,'\n%3s %5s ','P','coh');   fprintf(1,'%9.5f ',coeffs_coh_P);
fprintf(1,'\n%3s %5s ','P','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_P);
fprintf(1,'\n%3s %5s ','S','coh');   fprintf(1,'%9.5f ',coeffs_coh_S);
fprintf(1,'\n%3s %5s ','S','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_S);
fprintf(1,'\n%3s %5s ','Na','coh');   fprintf(1,'%9.5f ',coeffs_coh_Na);
fprintf(1,'\n%3s %5s ','Na','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_Na);
fprintf(1,'\n%3s %5s ','Cl','coh');   fprintf(1,'%9.5f ',coeffs_coh_Cl);
fprintf(1,'\n%3s %5s ','Cl','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_Cl);
fprintf(1,'\n%3s %5s ','Ar','coh');   fprintf(1,'%9.5f ',coeffs_coh_Ar);
fprintf(1,'\n%3s %5s ','Ar','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_Ar);
fprintf(1,'\n%3s %5s ','K','coh');   fprintf(1,'%9.5f ',coeffs_coh_K);
fprintf(1,'\n%3s %5s ','K','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_K);
fprintf(1,'\n%3s %5s ','Ca','coh');   fprintf(1,'%9.5f ',coeffs_coh_Ca);
fprintf(1,'\n%3s %5s ','Ca','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_Ca);
fprintf(1,'\n%3s %5s ','Fe','coh');   fprintf(1,'%9.5f ',coeffs_coh_Fe);
fprintf(1,'\n%3s %5s ','Fe','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_Fe);
fprintf(1,'\n%3s %5s ','Co','coh');   fprintf(1,'%9.5f ',coeffs_coh_Co);
fprintf(1,'\n%3s %5s ','Co','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_Co);
fprintf(1,'\n%3s %5s ','Cd','coh');   fprintf(1,'%9.5f ',coeffs_coh_Cd);
fprintf(1,'\n%3s %5s ','Cd','incoh'); fprintf(1,'%9.5f ',coeffs_incoh_Cd);
fprintf(1,'\n\n');