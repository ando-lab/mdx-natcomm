classdef Peptide
    % PEPTIDE - build amino acid chains from sequence
    %
    % development resources:
    %
    % * nomenclature: http://schmieder.fmp-berlin.info/teaching/educational_scripts/pdf/aminoacids.pdf
    % * amino acid chart: https://upload.wikimedia.org/wikipedia/commons/a/a9/Amino_Acids.svg
    % * chemical component dictionary: https://www.wwpdb.org/data/ccd/
    % * look up components in the dictionary: https://www.ebi.ac.uk/pdbe-srv/pdbechem/
    % * protparam: https://web.expasy.org/protparam/
    % * atomic weights: https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses
    properties
        sequence = ''; % e.g. 'EWRYKK'
        phi = []       % optional list of phi angles in degrees
        psi = []       % optional list of psi angles in degrees
        isCis = []     % optional logical list (empty means false: all angles trans)
        deprotonateHis = true; % if true, His will be neutral (otherwise positive)
        deprotonateCterm = true; % if true, C-term will be negatively charged
        deprotonateAsp = true; % if true, Asp will be negatively charged
        deprotonateGlu = true; % if true, Glu will be negatively charged
    end
    properties(Constant=true)
        peptideBondLength = 1.33; % constant (hard-coded at 1.33)
    end
    properties(Dependent=true)
        residues                  % (dependent) list of amino acids (mmCif IDs)
    end

    methods
        function obj = Peptide(seq)
            if nargin==0
                return;
            end

            if ~isempty(seq) && ischar(seq)
                obj.sequence = seq;

            elseif ~isempty(seq) && iscellstr(seq)
                obj.sequence = obj.aminolookup(seq);
            end

        end
        function val = get.residues(obj)
            res = obj.aminolookup(obj.sequence);
            nRes = length(res);

            % de-protonate C-terminus by default
            if obj.deprotonateCterm
                res{end} = [res{end} '_LEO2'];
            end

            % NH2 on N-terminus by defualt
            if nRes > 1
                res{1} = [res{1} '_LFOH'];
            end

            % use internal linking residues (type LL) by default
            if nRes > 2
                res(2:(end-1)) = cellfun(@(v) [v '_LL'],...
                    res(2:(end-1)),...
                    'UniformOutput',false);
            end

            if obj.deprotonateAsp
                % de-protonate ASP residues
                isASP = strncmpi(res,'asp',3);
                res(isASP) = cellfun(@(v) [v '_DHD2'],...
                    res(isASP),...
                    'UniformOutput',false);
            end

            if obj.deprotonateGlu
                % de-protonate GLU residues
                isGLU = strncmpi(res,'glu',3);
                res(isGLU) = cellfun(@(v) [v '_DHE2'],...
                    res(isGLU),...
                    'UniformOutput',false);
            end

            if obj.deprotonateHis
                % de-protonate His residues
                isHIS = strncmpi(res,'his',3);
                res(isHIS) = cellfun(@(v) [v '_DHD1'],...
                    res(isHIS),...
                    'UniformOutput',false);
            end
            % leave His neutral?
            val = res;
        end

        function M = build(obj)
            AminoAcids = model.chem.sequence2AminoAcids(obj.residues);
            for j=1:length(AminoAcids)
                if ~isempty(obj.phi)
                    AminoAcids(j).phi = obj.phi(j);
                end
                if ~isempty(obj.psi)
                    AminoAcids(j).psi = obj.psi(j);
                end
            end

            d = obj.peptideBondLength;

            PeptideBond = table();
            PeptideBond.id1 = {''};
            PeptideBond.id2 = {''}; % set later
            PeptideBond.isDouble = true; % has double-bond character (planar)
            PeptideBond.isAromatic = false;
            PeptideBond.isCovalent = true;

            for j=2:length(AminoAcids)
                M1 = AminoAcids(j-1);
                M2 = AminoAcids(j);

                zn = M1.c_peptide_normal;
                yn = M1.c_bond;
                xn = cross(yn,zn);
                R1 = [xn(:),yn(:),zn(:)];

                zn = M2.n_peptide_normal;
                yn = -M2.n_bond;
                xn = cross(yn,zn);

                if ~isempty(obj.isCis) && obj.isCis(j-1)
                    R2 = [xn(:),yn(:),zn(:)];
                else % trans
                    R2 = [-xn(:),yn(:),-zn(:)];
                end

                M2 = M2.translate(-M2.n_term - d*M2.n_bond).rotate(R2');
                M2 = M2.rotate(R1).translate(M1.c_term);
                M2 = M2.setCharge('N',0);
                M1 = M1.setCharge('C',0);
                if ~M2.isProline
                    AminoAcids(j) = M2.removeAtom('HN2');
                else
                    AminoAcids(j) = M2.removeAtom('H2');
                end
                AminoAcids(j-1) = M1.removeAtom({'HXT','OXT'});
            end
            for j=1:length(AminoAcids)
                AminoAcids(j) = AminoAcids(j).prependToName([num2str(j) '/']);
            end
            M = model.chem.Molecule();
            M = M.addMolecule(AminoAcids(1));
            for j=2:length(AminoAcids)
                M = M.addMolecule(AminoAcids(j));
                PeptideBond.id1 = [num2str(j-1),'/C'];
                PeptideBond.id2 = [num2str(j),'/N'];
                M = M.addBond(PeptideBond);
            end
            M.name = {AminoAcids.name};
            M.id = {AminoAcids.id};
        end
    end
    methods(Static)
        function aa = aminolookup(seq)
            % Similar functionality to MATLAB's aminolookup function for
            % translating between single-letter and three-letter
            % abbreviations. More complex lookups are not supported.
            %
            % The lookup tables were generated automatically using MATLAB's
            % aminolookup function (run the following code)
            %
            %     l1 = 'ACDEFGHIKLMNPQRSTVWY';
            %     fprintf(1,'l3lookup = struct(');
            %     for j=1:length(l1)
            %         fprintf(1,'''%s'',''%s'',',l1(j),aminolookup(l1(j)));
            %     end
            %     fprintf(1,'\b);\n');
            %     fprintf(1,'l1lookup = struct(');
            %     for j=1:length(l1)
            %         fprintf(1,'''%s'',''%s'',',lower(aminolookup(l1(j))),l1(j));
            %     end
            %     fprintf(1,'\b);\n');

            l3lookup = struct('A','Ala','C','Cys','D','Asp','E','Glu','F','Phe','G','Gly','H','His','I','Ile','K','Lys','L','Leu','M','Met','N','Asn','P','Pro','Q','Gln','R','Arg','S','Ser','T','Thr','V','Val','W','Trp','Y','Tyr');
            l1lookup = struct('ala','A','cys','C','asp','D','glu','E','phe','F','gly','G','his','H','ile','I','lys','K','leu','L','met','M','asn','N','pro','P','gln','Q','arg','R','ser','S','thr','T','val','V','trp','W','tyr','Y');

            if ischar(seq)
                aa = cell(1,length(seq));
                for j=1:length(aa)
                    aa{j} = l3lookup.(upper(seq(j)));
                end
            elseif iscellstr(seq)
                for j=1:length(seq)
                    seq{j} = l1lookup.(lower(seq{j}));
                end
                aa = strjoin(seq,'');
            else
                aa = {};
            end
        end
    end
end
