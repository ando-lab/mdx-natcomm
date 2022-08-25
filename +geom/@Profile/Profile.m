classdef Profile
    %PROFILE a collection of methods for spot profile calculations
    %   <Detailed explanation goes here>
    
    methods(Static)
        
        function [idx,hkl_refl] = hkl2idx(h,k,l)
            % assign each pixel an index to its nearest reflection
            h = round(h);
            k = round(k);
            l = round(l);
            [hkl_refl,~,ic] = unique([h(:),k(:),l(:)],'rows');
            idx = reshape(ic,size(h));
        end
        
        function [e1,e2,e3] = s2e(s_refl,Source)
            % calculate local coordinate system for each reflection
            s0 = Source.wavevector;
            
            e1 = zeros(size(s_refl));
            e2 = zeros(size(s_refl));
            e3 = zeros(size(s_refl));
            
            for j=1:size(s_refl,1)
                k1 = s0 + s_refl(j,:);
                e1(j,:) = quickcross(k1,s0);
                e1(j,:) = e1(j,:)/sqrt(quickdot(e1(j,:),e1(j,:)));
                
                e2(j,:) = quickcross(k1,e1(j,:));
                e2(j,:) = e2(j,:)/sqrt(quickdot(e2(j,:),e2(j,:)));
                
                e3(j,:) = (k1 + s0);
                e3(j,:) = e3(j,:)/sqrt(quickdot(e3(j,:),e3(j,:)));
            end
        end
        
        function [phixy1,phixy2] = s2phi(s_refl,Source,Spindle,Detector)
            % I changed how this works on 3/25/2019. It used to return:
            % [phi, x, y] for whichever solution had the lowest value of
            % |phi|.
            
            s0 = Source.wavevector;
            m2 = Spindle.rotationAxis;
            f = Detector.f;
            ed = Detector.ed;
            [phixy1, phixy2] = calc_reflection_center(s0, m2, f, ed, s_refl);
            
%            [c1,c2] = calc_reflection_center(s0, m2, f, ed, s_refl);
%             all_phi = [c1(:,1),c2(:,1)];
%             [~,ix] = min(abs(all_phi),[],2);
%             
%             c = c1;
%             for j=1:size(all_phi,1)
%                 if ix(j)==2
%                     c(j,:) = c2(j,:);
%                 end
%             end
%             
%             x = c(:,2);
%             y = c(:,3);
%             phi = c(:,1);
        end
        
        function [eps1,eps2] = s2eps(sx,sy,sz,idx,s,e1,e2,Source,Detector)
            
            nrefl = size(s,1);
            
            s0 = Source.wavevector;
            ed = Detector.ed;
            
            d1 = ed(:,1)';
            d2 = ed(:,2)';
            d3 = ed(:,3)';
            
            krefl = s + repmat(s0,nrefl,1);
            krefl_mag = sqrt(sum(krefl.*krefl,2));
            
            d1_dot_e1 = sum(repmat(d1,nrefl,1).*e1,2);
            d1_dot_e2 = sum(repmat(d1,nrefl,1).*e2,2);
            
            d2_dot_e1 = sum(repmat(d2,nrefl,1).*e1,2);
            d2_dot_e2 = sum(repmat(d2,nrefl,1).*e2,2);
            
            d3_dot_e1 = sum(repmat(d3,nrefl,1).*e1,2);
            d3_dot_e2 = sum(repmat(d3,nrefl,1).*e2,2);

            sx_refl = s(:,1);
            sy_refl = s(:,2);
            sz_refl = s(:,3);
            
            eps1 = (180/pi)*((sx - sx_refl(idx)).*d1_dot_e1(idx) + ...
                (sy - sy_refl(idx)).*d2_dot_e1(idx) + ...
                (sz - sz_refl(idx)).*d3_dot_e1(idx))./krefl_mag(idx);
            
            eps2 =(180/pi)* ((sx - sx_refl(idx)).*d1_dot_e2(idx) + ...
                (sy - sy_refl(idx)).*d2_dot_e2(idx) + ...
                (sz - sz_refl(idx)).*d3_dot_e2(idx))./krefl_mag(idx);

        end
        
        function Rj = partiality(phi,e1,sigmaM,Spindle,frame)
            % added frame argument on March 25, 2019
            m2 = Spindle.rotationAxis;
            dPhi = Spindle.oscillationRange;
            phiMid = Spindle.frame2phi(frame);
            phi1 = phiMid - dPhi/2;
            phi2 = phiMid + dPhi/2; % before, assumed phimid = 0
            %disp([phi1,phi2])
            Rj = calc_partiality(phi,e1,phi1,phi2,sigmaM,m2);
           
        end
        
    end
end

function c = quickdot(a,b)
c = a(1)*b(1) + a(2)*b(2) + a(3)*b(3);
end

function c = quickcross(a,b)
c = zeros(1,3);
c(1) = a(2)*b(3) - a(3)*b(2);
c(2) = a(3)*b(1) - a(1)*b(3);
c(3) = a(1)*b(2) - a(2)*b(1);
end

