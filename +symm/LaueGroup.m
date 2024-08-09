classdef LaueGroup < symm.PointGroup
    % LaueGroup - symmetry operators and other information related to a
    % crystallographic Laue group
    %
    % Implemented as a set of dynamic properties that wraps
    % around a SpaceGroupInfo object
    %
    % Also includes the testASU property, which is a function handle that
    % returns true for all h,k,l within the reference asu.
    
    properties(Dependent=true)
        testASU
    end
    
    methods
        function obj = LaueGroup(varargin)
            %LaueGroup Construct an instance of this class
            obj = obj@symm.PointGroup(varargin{:});
            assert(obj.Info.isLaueGroup,...
                'requested space group is not a Laue group');
        end
        
        function isasu = get.testASU(obj)
            % TESTASU - Test whether a given miller index is in the reference asymmetric unit
            %
            % The definitions are based on the conventions used in cctbx
            % See: cctbx/sgtbx/reciprocal_space_ref_asu.cpp
            switch obj.Info.number
                case 2 % -1
                    isasu = @(h,k,l) l>0 | (l==0 & (h>0 | (h==0 & k>=0)));
                case 10 % 2/m
                    isasu = @(h,k,l) k>=0 & (l>0 | (l==0 & h>=0));
                case 47 % mmm
                    isasu = @(h,k,l) h>=0 & k>=0 & l>=0;
                case {83,175} % 4/m or 6/m
                    isasu = @(h,k,l) l>=0 & ((h>=0 & k>0) | (h==0 & k==0));
                case {123,191} % 4/mmm or 6/mmm
                    isasu = @(h,k,l) h>=k & k>=0 & l>=0;
                case 147 % -3
                    isasu = @(h,k,l) (h>=0 & k>0) | (h==0 & k==0 & l>=0);
                case 162 % -31m
                    isasu = @(h,k,l) h>=k & k>=0 & (k>0 | l>=0);
                case 164 % -3m1
                    isasu = @(h,k,l) h>=k & k>=0 & (h>k | l>=0);
                case 166 % -3m(R)
                    isasu = @(h,k,l) h>=0 & k>=h & l<=h & ~(h==0 & l<0 & abs(l)>k);
                case 200 % m-3
                    isasu = @(h,k,l) h>=0 & ((l>=h & k>h) | (l==h & k==h));
                case 221 % m-3m
                    isasu = @(h,k,l) k>=l & l>=h & h>=0;
                otherwise
                    error('Laue group not found... this should not happen!');
            end
        end
    end
end

%


