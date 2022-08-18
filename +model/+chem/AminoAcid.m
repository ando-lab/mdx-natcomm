classdef AminoAcid < model.chem.Molecule
    %AMINOACID Coordinates of an amino acid model
    
    properties(Dependent = true)
        n_term
        c_term
        n_peptide_normal
        c_peptide_normal
        ca_normal
        n_bond
        c_bond
        n_term_atoms
        c_term_atoms
        phi
        psi
        isProline
    end
    
    methods
        function obj = AminoAcid(varargin)
            
            if ~isempty(varargin) && isa(varargin{1},'model.chem.Molecule')
                args = {};
                isMolecule = true;
            else
                args = varargin;
                isMolecule = false;
            end
            obj@model.chem.Molecule(args{:});
            
            if isMolecule
                obj.G = varargin{1}.G;
                obj.name = varargin{1}.name;
                obj.id = varargin{1}.id;
            end
        end
        function val = get.isProline(obj)
            val = strncmpi(obj.id,'pro',3);
        end
        function val = get.c_term_atoms(obj)
            G0 = obj.G.rmedge('C','CA');
            val = G0.bfsearch('C');
        end
        function val = get.n_term_atoms(obj)
            if ~obj.isProline
                G0 = obj.G.rmedge('N','CA');
                val = G0.bfsearch('N');
            else % it's proline... no changes allowed
                val = {};
            end
        end
        
        function val = get.c_term(obj)
            c = obj.findAtom('C');
            val = [c.x,c.y,c.z];
        end
        function val = get.n_term(obj)
            n = obj.findAtom('N');
            val = [n.x,n.y,n.z];
        end
        function val = get.ca_normal(obj)
            c = obj.findAtom('C');
            ca = obj.findAtom('CA');
            n = obj.findAtom('N');
            c_coords = [c.x,c.y,c.z];
            ca_coords = [ca.x,ca.y,ca.z];
            n_coords = [n.x,n.y,n.z];
            can = cross(c_coords-ca_coords,ca_coords-n_coords);
            val = can/sqrt(dot(can,can));
        end
        function val = get.c_peptide_normal(obj)
            c = obj.findAtom('C');
            ca = obj.findAtom('CA');
            o = obj.findAtom('O');
            c_coords = [c.x,c.y,c.z];
            ca_coords = [ca.x,ca.y,ca.z];
            o_coords = [o.x,o.y,o.z];
            cpn = cross(ca_coords-o_coords,c_coords-ca_coords);
            val = cpn/sqrt(dot(cpn,cpn));
        end
        
        function val = get.n_peptide_normal(obj)
            if ~obj.isProline
                n = obj.findAtom('N');
                h = obj.findAtom('H');
                ca = obj.findAtom('CA');
                
                n_coords = [n.x,n.y,n.z];
                h_coords = [h.x,h.y,h.z];
                ca_coords = [ca.x,ca.y,ca.z];
                
                npn = -1*cross(ca_coords-n_coords,n_coords-h_coords);
                val = npn/sqrt(dot(npn,npn));
            else % it's proline... peptide normal fixed so that phi is -60
                angle = 60*pi/180;
                
                n = obj.findAtom('N');
                ca = obj.findAtom('CA');
                can = obj.ca_normal;
                rot_axis = [ca.x - n.x,ca.y - n.y,ca.z - n.z];
                zn = rot_axis/sqrt(dot(rot_axis,rot_axis));
                
                xn = can;
                yn = cross(zn,xn);
                R1 = [xn(:),yn(:),zn(:)];
                
                R2 = [cos(angle),-sin(angle),0; sin(angle),cos(angle),0; 0,0,1];
                
                val = (R1*R2*R1'*can(:))';
            end
        end
        
        function val = get.c_bond(obj)
            % angle between C-CA and C-N bonds
            angle_ca = 116*pi/180;
            
            c = obj.findAtom('C');
            c_coords = [c.x,c.y,c.z];
            
            ca = obj.findAtom('CA');
            ca_coords = [ca.x,ca.y,ca.z];
            
            pcn = obj.c_peptide_normal;
            
            xn = ca_coords - c_coords;
            xn = xn/sqrt(dot(xn,xn));
            
            yn = cross(pcn,xn);
            
            val = xn*cos(angle_ca) + yn*sin(angle_ca);
        end
        
        function val = get.n_bond(obj)
            % angle between N-CA and N-C bonds
            angle_ca = 122*pi/180;
            
            n = obj.findAtom('N');
            n_coords = [n.x,n.y,n.z];
            
            ca = obj.findAtom('CA');
            ca_coords = [ca.x,ca.y,ca.z];
            
            pnn = obj.n_peptide_normal;
            
            xn = ca_coords - n_coords;
            xn = xn/sqrt(dot(xn,xn));
            
            yn = -cross(pnn,xn);
            
            val = xn*cos(angle_ca) + yn*sin(angle_ca);
        end
        
        function val = get.psi(obj)
            pcn = obj.c_peptide_normal;
            can = obj.ca_normal;
            %val = (180/pi)*acos(dot(pcn,can));
            
            ca = obj.findAtom('CA');
            c = obj.findAtom('C');
            bond_vec = [c.x-ca.x,c.y-ca.y,c.z-ca.z];
            bond_vec = bond_vec/sqrt(dot(bond_vec,bond_vec));
            val = (180/pi)*atan2(dot(bond_vec,cross(can,pcn)),dot(pcn,can));
        end
        
        function obj = set.psi(obj,psiNew)
            
            angle = pi*(psiNew - obj.psi)/180;
            
            term_idx = obj.G.findnode(obj.c_term_atoms);
            all_coords = table2array(obj.G.Nodes(term_idx,{'x','y','z'}))';
            
            c = obj.findAtom('C');
            ca = obj.findAtom('CA');
            
            rot_axis = [c.x-ca.x,c.y-ca.y,c.z-ca.z];  % c -> ca
            rot_axis = rot_axis/sqrt(dot(rot_axis,rot_axis));
            
            zn = rot_axis;
            xn = obj.c_peptide_normal;
            yn = cross(zn,xn);
            
            R1 = [xn(:),yn(:),zn(:)];
            
            R2 = [cos(angle),-sin(angle),0; sin(angle),cos(angle),0; 0,0,1];
            
            c1 = R1*R2*R1'*(all_coords - repmat([c.x;c.y;c.z],1,size(all_coords,2)));
            
            obj.G.Nodes.x(term_idx) = c1(1,:)' + c.x;
            obj.G.Nodes.y(term_idx) = c1(2,:)' + c.y;
            obj.G.Nodes.z(term_idx) = c1(3,:)' + c.z;
            
        end
        
        function val = get.phi(obj)
            pnn = obj.n_peptide_normal;
            can = obj.ca_normal;
            ca = obj.findAtom('CA');
            n = obj.findAtom('N');
            bond_vec = [ca.x-n.x,ca.y-n.y,ca.z-n.z];
            bond_vec = bond_vec/sqrt(dot(bond_vec,bond_vec));
            val = (180/pi)*atan2(dot(bond_vec,cross(pnn,can)),dot(pnn,can));
        end
        
        function obj = set.phi(obj,phiNew)
            if obj.isProline
                warning('cannot set phi angle of Proline');
                return;
            end
            
            angle = -pi*(phiNew - obj.phi)/180;
            
            term_idx = obj.G.findnode(obj.n_term_atoms);
            all_coords = table2array(obj.G.Nodes(term_idx,{'x','y','z'}))';
            
            ca = obj.findAtom('CA');
            n = obj.findAtom('N');
            
            rot_axis = [ca.x-n.x,ca.y-n.y,ca.z-n.z];  % ca -> n
            rot_axis = rot_axis/sqrt(dot(rot_axis,rot_axis));
            
            zn = rot_axis;
            xn = obj.n_peptide_normal;
            yn = cross(zn,xn);
            
            R1 = [xn(:),yn(:),zn(:)];
            
            % clockwise rotation (right handed)
            R2 = [cos(angle),-sin(angle),0; sin(angle),cos(angle),0; 0,0,1];
            
            c1 = R1*R2*R1'*(all_coords - repmat([n.x;n.y;n.z],1,size(all_coords,2)));
            
            obj.G.Nodes.x(term_idx) = c1(1,:)' + n.x;
            obj.G.Nodes.y(term_idx) = c1(2,:)' + n.y;
            obj.G.Nodes.z(term_idx) = c1(3,:)' + n.z;
            
        end
        
        function plot(obj)
            plot@model.chem.Molecule(obj);
            hold on;
            pcn = obj.c_peptide_normal;
            pnn = obj.n_peptide_normal;
            cbond = obj.c_bond;
            nbond = obj.n_bond;
            c = obj.c_term;
            n = obj.n_term;
            yy = [c(:),c(:) + pcn(:),c(:),c(:) + cbond(:)];
            plot3(yy(1,:),yy(2,:),yy(3,:),'c-','LineWidth',3);
            yy = [n(:),n(:) + pnn(:),n(:),n(:) + nbond(:)];
            plot3(yy(1,:),yy(2,:),yy(3,:),'r-','LineWidth',3);
        end
    end
end

