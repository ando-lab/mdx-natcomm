classdef DeltaPDFTools < util.propertyValueConstructor
    %DeltaPDFTools - 
    
    properties
        MT_peak % proc.script.MapTools object describing the full Patterson map
        rmax % maximum radius for peak extraction
        supercell
        smin = 0.2; % 1/d range for fitting in reciprocal space
        smax = 0.6;
        verbose = true;
    end
    properties(Dependent)
        G_sup
        peak_mask
    end
    
    methods
        function obj = DeltaPDFTools(varargin)
            %DeltaPDFTools
            obj@util.propertyValueConstructor(varargin{:});
        end
        
        function G_sup = get.G_sup(obj)
            G_sup = latt.PeriodicGrid(obj.supercell,[0,0,0],[1,1,1]).invert;
        end
        
        function peak_mask = get.peak_mask(obj)
            peak_mask = obj.MT_peak.spherical_mask(obj.rmax);
        end
        
        function d = read_peak_data(obj,h5In,mapname,dset)
            [n1,n2,n3] = obj.G_sup.grid();
            msk = obj.peak_mask;
            
            d = cell(obj.supercell);
            
            for j=1:prod(obj.supercell)
                if obj.verbose
                    fprintf('reading peak %d of %d\n',j,prod(obj.supercell));
                end
                MTj = obj.MT_peak;
                MTj.Grid.ori = MTj.Grid.ori + [n1(j),n2(j),n3(j)];
                d{j} = MTj.read_data(h5In,mapname,dset);
                d{j}(~msk) = NaN;
            end
            
        end
        
        function pdf_calc = simulate_patterson_peak(obj,AtomFF)
            
            MT_calc = obj.MT_peak;
            MT_calc.Grid = resize2primitive(obj.MT_peak.Grid);
            
            Fatom =proc.script.StructureFactor(...
                'LatticeGrid',latt.LatticeGrid(MT_calc.Grid,MT_calc.Basis),...
                'SpaceGroup',MT_calc.SpaceGroup,...
                'Badd',5, ... % <-- might need to play with this. 5 seems OK
                'rmax',4 ... % <-- 3 is default. 4 is conservative but kinda slow
                ).run(AtomFF);
            
            Icalc = real(Fatom.*conj(Fatom));
            
            [pdf_calc,MT_calc] = MT_calc.invert.fourier_transform(Icalc);
            
            [~,resizefun] = MT_calc.resize('grid',obj.MT_peak.Grid);
            pdf_calc = resizefun(pdf_calc);
            pdf_calc(~obj.peak_mask) = NaN;
        end
        
        
        function [fit_info,dpdf_fit] = deconvolve_peaks(obj,peak_data,peak_sim)
            
            [param2struct,v0] = obj.parameterize();
            
            [Ip_sim,MTr] = obj.MT_peak.fourier_transform(peak_sim);
            Ip_sim = double(Ip_sim);
            isIncl = MTr.spherical_mask([obj.smin,obj.smax]);
            
            fit_info = struct('V',[],'cc',[],'p',[]);
            fit_info = repmat(fit_info,size(peak_data));
            
            [f1,f2,f3] = MTr.Grid.grid();
            [sx,sy,sz] = MTr.Basis.frac2lab(f1,f2,f3);
            
            q = 2*pi*[sx(isIncl),sy(isIncl),sz(isIncl)]';
            
            Wfun = @(p) dot(q,p.V*q);
            fitfun = @(v,Iref,Imeas) (Iref.*reshape(Wfun(param2struct(v)),size(Iref)) - Imeas);
            
            for j=1:numel(peak_data)
                if obj.verbose
                    fprintf(1,'fitting peak %d of %d\n',j,numel(peak_data));
                end
                Ip_meas = obj.MT_peak.fourier_transform(peak_data{j});
                Ip_meas = double(Ip_meas);
                xfit = lsqnonlin(@(v) fitfun(v,Ip_sim(isIncl),Ip_meas(isIncl)),v0);
                
                % calculate cc
                xfitparam = param2struct(xfit);
                Ip_fit = Ip_sim(isIncl).*Wfun(xfitparam)';
                cc = sum(Ip_fit.*Ip_meas(isIncl))./sqrt(sum(Ip_fit.^2)*sum(Ip_meas(isIncl).^2));
                
                % save outputs
                fit_info(j).p = xfit;
                fit_info(j).V = xfitparam.V;
                fit_info(j).cc = cc;
                
            end
            
            
            % compute the fitted peaks
            dpdf_fit = cell(size(peak_data));
            
            q = 2*pi*[sx(:),sy(:),sz(:)]';
            
            for j=1:numel(peak_data)
                wfit = reshape(dot(q,fit_info(j).V*q),size(sx));
                dpdf_fit{j} = MTr.fourier_transform(Ip_sim.*wfit);
            end
            
            
        end
        
        
        
        function [p2s,v0] = parameterize(obj)
            
            % note... I should think more carefully about the symmetries of
            % V. But I think for now it's OK to assume that the data will
            % tell me the symmetry. I can just use P1.
            
            symm = 'triclinic';
            %symm = obj.MT_peak.SpaceGroup.Info.latticeSystem;
            
            switch lower(symm)
                case 'cubic'
                    p2s = @(v) struct(...
                        'V',[v(1),0,0;...
                             0,v(1),0;...
                             0,0,v(1)]);
                    v0 = [.1];
                case {'hexagonal','rhombohedral'}
                    
                    p2s = @(v) struct(...
                        'V',[v(1),v(1),0;...
                             v(1),v(1),0;...
                             0,0,v(2)]);
                    v0 = [.1,.1];
                    
                case 'monoclinic'
                    
                    p2s = @(v) struct(...
                        'V',[v(1),0,v(4);...
                             0,v(2),0;...
                             v(4),0,v(3)]);
                    v0 = [.1,.1,.1,0];
                    
                case 'orthorhombic'
                    p2s = @(v) struct(...
                        'V',[v(1),0,0;...
                             0,v(2),0;...
                             0,0,v(3)]);
                    v0 = [.1,.1,.1];
                case 'tetragonal'
                    p2s = @(v) struct(...
                        'V',[v(1),0,0;...
                             0,v(1),0;...
                             0,0,v(2)]);
                    v0 = [.1,.1];
                case 'triclinic'
                    p2s = @(v) struct(...
                        'V',[v(1),v(4),v(5);...
                             v(4),v(2),v(6);...
                             v(5),v(6),v(3)]);
                    v0 = [.1,.1,.1,0,0,0];
                otherwise
                    error('oops');
            end
            
            
        end
    end
    
    
    
    methods(Static)

    end
end






function G0 = resize2primitive(G)

G0 = G;
G0.N = round(G0.invert.P);
G0.P = [1,1,1];
G0 = G0.invert.invert; % recenter
end