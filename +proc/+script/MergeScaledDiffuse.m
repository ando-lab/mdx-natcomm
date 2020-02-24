classdef MergeScaledDiffuse < util.propertyValueConstructor
    properties
        Grid = grid.Sub3d('ndiv',[1,1,1]);
        Crystal = geom.Crystal.empty();
        tmpDir = ''; % gets set in constructor
        makeTmpDir = true; % will make and delete the temporary directory if true
        workingDirectory = './'
        fid = 1; % first argument of fprintf (1 = comannd window)
        nMin = 3; % minimum multiplicity for merging (n < nMin are discarded)
        residMax = 4; % observations with residual > residMax*sigma are discarded
    end
    
    methods
        function obj = MergeScaledDiffuse(varargin)
            obj@util.propertyValueConstructor(varargin{:});
            if isempty(obj.tmpDir) && obj.makeTmpDir
                obj.tmpDir = fullfile(tempdir,...
                    sprintf('mdx-kit-%s',datestr(now,'yyyymmddHHMMSS')));
            end
        end
        
        function fn = mapToColumns(obj,inputFileNames)

            % set up
            fprintf(obj.fid,'START\n');
            
            % set temporary file names
            if ~isdir(obj.tmpDir)
                if obj.makeTmpDir
                    fprintf(obj.fid,'making temporary directory %s\n',obj.tmpDir);
                    mkdir(obj.tmpDir);
                else
                    error('temporary directory does not exist');
                end
            else
                if obj.makeTmpDir
                    error('temporary directory already exists but makeTmpDir = true.');
                else
                    fprintf(obj.fid,'using temporary directory %s\n',obj.tmpDir);
                end
            end
            
            ncols = prod(obj.Grid.ndiv);
            fn = cell(1,ncols);
            % set names for temporary files (one for each column)
            for j=1:ncols
                fn{j} = fullfile(obj.tmpDir,sprintf('%d.dat',j));
            end
            
            % read in diffuseTable files and write temporary files
            for dirIndex=1:length(inputFileNames)
                thisInputFileName = fullfile(obj.workingDirectory,inputFileNames{dirIndex});
                fprintf(obj.fid,'loading diffuseTable from %s\n',thisInputFileName);
                % load diffuse table
                inputData = load(thisInputFileName,'diffuseTable');
                
                fprintf(obj.fid,'  mapping hkl to asu\n');
                % map hkl to asu
                [h,k,l] = obj.Grid.hkl2hkl(...
                    inputData.diffuseTable.h,...
                    inputData.diffuseTable.k,...
                    inputData.diffuseTable.l,...
                    inputData.diffuseTable.dh,...
                    inputData.diffuseTable.dk,...
                    inputData.diffuseTable.dl);
                [h,k,l] = obj.Crystal.hkl2asu(h,k,l);
                [sub1,sub2] = obj.Grid.hkl2index(h,k,l);
                clear h k l % make space
                
                fprintf(obj.fid,'  get rows of diffuseTable to put in each binary file\n');
                [~,ix] = ismember(sub2,1:ncols);
                ib = accumarray(ix,1:numel(ix),[ncols,1],@(x) {x});
                clear ix sub2 % make space
                
                fprintf(obj.fid,'  append to binary files\n');
                for j=1:ncols
                    A = [sub1(ib{j}),...
                        inputData.diffuseTable.I(ib{j}),...
                        inputData.diffuseTable.sigma(ib{j})];
                    thisFileID = fopen(fn{j},'a');
                    fwrite(thisFileID,A','double'); % transposed A because fread can only have indeterminant cols
                    fclose(thisFileID);
                end
                clear A diffuseTable sub1 ib
                
            end
            
            fd = cellfun(@dir,fn);
            totalMb = sum([fd.bytes])/1E6;
            fprintf(obj.fid,'DONE: %d tmp files written (%.3f Mb total)\n',length(fd),totalMb); % 2355.343 Mb written
            
        end
        
        function hklTable = mergeRandomHalfSets(obj,fn,isincl)
            % proc.script.mergeRandomHalfSets
            %
            % random split is done by proc.scale.weightedRandomSplit
            %
            % note: does not do outlier detection. run mergeColumns first, and
            % pass isincl as an argument.
            fprintf(obj.fid,'START: mergeRandomHalfSets\n\n');
            
            ncols = prod(obj.Grid.ndiv);
            varNames = {'h','k','l','dh','dk','dl','I1','sigma1','I2','sigma2'};
            allTables = cell(ncols,1);
            
            for j=1:ncols
                fprintf(obj.fid,'processing file %d of %d\n',j,ncols);
                
                % read data from temporary file
                thisFileID = fopen(fn{j},'r');
                A = fread(thisFileID,[3,Inf],'double')';
                fclose(thisFileID);
                if isempty(A) % nothing in this voxel group
                    allTables{j} = zeros(0,length(varNames));
                    %isIn1{j} = false(0,1);
                    continue;
                end
                
                sub1 = A(:,1);
                I = A(:,2);
                sigma = A(:,3);
                
                if nargin > 2 && ~isempty(isincl)
                    sub1 = sub1(isincl{j});
                    I = I(isincl{j});
                    sigma = sigma(isincl{j});
                end
                
                % merge reflections
                [subRef,~,ih] = unique(sub1);
                hMax = length(subRef);
                
                isgrp1 = proc.scale.weightedRandomSplit(ih,1./sigma.^2);
                
                [Imerge1,sigmaMerge1] = mergeAndReject(...
                    ih(isgrp1),hMax,I(isgrp1),sigma(isgrp1));
                [Imerge2,sigmaMerge2] = mergeAndReject(...
                    ih(~isgrp1),hMax,I(~isgrp1),sigma(~isgrp1));
                
                % make output table
                [h,k,l] = obj.Grid.index2hkl(subRef);
                [~,~,~,dh,dk,dl] = obj.Grid.index2hkl(1,j);
                dh = repmat(dh,size(h,1),1);
                dk = repmat(dk,size(h,1),1);
                dl = repmat(dl,size(h,1),1);
                
                allTables{j} = [h,k,l,dh,dk,dl,Imerge1,sigmaMerge1,Imerge2,sigmaMerge2];
                
            end
            fprintf(obj.fid,'constructing hklTable\n');
            hklTable = array2table(cell2mat(allTables),'VariableNames',varNames);
            clear allTables
            fprintf(obj.fid,'sorting hklTable\n');
            hklTable = sortrows(hklTable,{'h','k','l','dh','dk','dl'});
            
            fprintf(obj.fid,'DONE\n');
        end
        
        function [hklTable,isincl] = mergeColumns(obj,fn)
            fprintf(obj.fid,'START: mergeColumns\n\n');
            
            ncols = prod(obj.Grid.ndiv);
            varNames = {'h','k','l','dh','dk','dl','I','sigma','n','x2'};
            allTables = cell(ncols,1);
            isincl = cell(ncols,1);
            
            ntot = 0;
            nreject = 0;
            
            % read in the temporary files
            for j=1:ncols
                fprintf(obj.fid,'processing file %d of %d\n',j,ncols);
                
                thisFileID = fopen(fn{j},'r');
                A = fread(thisFileID,[3,Inf],'double')';
                fclose(thisFileID);
                if isempty(A) % nothing in this voxel group
                    allTables{j} = zeros(0,length(varNames));
                    isincl{j} = false(0,1);
                    continue;
                end
                
                sub1 = A(:,1);
                I = A(:,2);
                sigma = A(:,3);
                %clear A
                
                % merge reflections
                [subRef,~,ih] = unique(sub1);
                hMax = length(subRef);
                
                [Imerge,sigmaMerge,n,x2,isincl{j}] = mergeAndReject(...
                    ih,hMax,I,sigma,obj.residMax,obj.nMin);
                
                nreject = nreject + sum(~isincl{j});
                ntot = ntot + size(A,1);
                
                % make output table
                [h,k,l] = obj.Grid.index2hkl(subRef);
                [~,~,~,dh,dk,dl] = obj.Grid.index2hkl(1,j);
                dh = repmat(dh,size(h,1),1);
                dk = repmat(dk,size(h,1),1);
                dl = repmat(dl,size(h,1),1);
                
                allTables{j} = [h,k,l,dh,dk,dl,Imerge,sigmaMerge,n,x2];
                %hklMerge = [hklMerge; t];
            end %
            fprintf(obj.fid,'constructing hklTable\n');
            hklTable = array2table(cell2mat(allTables),'VariableNames',varNames);
            clear allTables
            fprintf(obj.fid,'sorting hklTable\n');
            hklTable = sortrows(hklTable,{'h','k','l','dh','dk','dl'});
            fprintf(obj.fid,'DONE: %d reflections, %d rejected, merged to %d\n',ntot,nreject,size(hklTable,1));
        end
        
        function tf = clearTmp(obj,fn)
            
            % delete temporary files and directory
            cellfun(@delete,fn);
            
            % remove temporary directory if it was created
            if obj.makeTmpDir
                rmdir(obj.tmpDir);
            end
            
            tf = true;
        end 
    end
end


function [Imerge,sigmaMerge,n,x2,isincl] = mergeAndReject(...
    ih,hMax,I,sigma,residMax,nMin)

if nargin<5 || isempty(residMax)
    residMax = Inf;
end
if nargin<6 || isempty(nMin)
    nMin = 0;
end

w0 = 1./sigma.^2;
w = accumarray(ih,w0,[hMax,1]);
Imerge = accumarray(ih,I.*w0,[hMax,1])./w;
sigmaMerge = sqrt(1./w);

% check outliers
resid = (Imerge(ih) - I)./sigma;
n = accumarray(ih(~isinf(sigma)),1,[hMax,1]); % how many times each voxel was measured

if ~isinf(residMax) || nMin > 0 % need to apply filter
    % apply merge filtering rules:
    % - a voxel must be measured nMin or more times initially
    % - the absolute value of the residual must be < residMax*sigma
    sigma(abs(resid)>residMax | n(ih) < nMin) = Inf;
    
    % merge again
    w0 = 1./sigma.^2;
    w = accumarray(ih,w0,[hMax,1]);
    Imerge = accumarray(ih,I.*w0,[hMax,1])./w;
    sigmaMerge = sqrt(1./w);
    
    % calculate chi-squared
    resid = (Imerge(ih) - I)./sigma;
    n = accumarray(ih(~isinf(sigma)),1,[hMax,1]); % how many times each voxel was measured
end

isincl = ~isinf(sigma);
x2 = accumarray(ih,resid.^2,[hMax,1]);
end


