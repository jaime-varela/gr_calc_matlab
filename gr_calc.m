classdef gr_calc  < handle
    %gr_calc a class to compute standard GR objects
    %   constructor: gr_calc(metric,coordinates)
    %   various methods to compute standard general relativity objects.
    
    properties
        grMetric;
        grCoordinates;
        grDimension;
        grIMetric;
        grMetricDerv = {};
        grIMetricDerv = {};

        % usage \Gamma^i_{kl} = this.grChristoffel(k,l,i)
        grChristoffel = {};
        grRiemann = {};
        grRicci = {};
        grRscalar = {};
        grEinstein = {};
    end

    properties (Access = private)
        % bool variables representing wether or not 
        % some gr objects have been computed
        isComputedChristoffel = false;
        isComputedRiemann = false;
        isComputedRicci = false;
        isComputedRScalar = false;
        isComputedEinstein = false;
    end
    
    methods
        function obj = gr_calc(symbolicMetricMatrix,symbCoordinateVector)
            %gr_calc Load a symbolic metric and its coordinates
            %   Indexes start at one, keeping with Matlab's convention
            obj.grMetric = symbolicMetricMatrix;
            obj.grCoordinates = symbCoordinateVector;
            obj.grDimension = length(symbCoordinateVector);
            % calculated quantities
            obj.grIMetric = inv(symbolicMetricMatrix);
            % derivatives of the metric
            obj.grMetricDerv = sym('grStuff',[obj.grDimension,obj.grDimension,obj.grDimension]);
            for coordDerv = 1:length(symbCoordinateVector)
                obj.grMetricDerv(:,:,coordDerv) = diff(obj.grMetric,symbCoordinateVector(coordDerv));
            end
            obj.grIMetricDerv = sym('grStuff',[obj.grDimension,obj.grDimension,obj.grDimension]);
            for iCoordDerv = 1:length(symbCoordinateVector)
               obj.grIMetricDerv(:,:,iCoordDerv) = -1*(obj.grIMetric)*( obj.grMetricDerv(:,:,iCoordDerv) )*(obj.grIMetric);
            end
        end
        
        function christoffelIndexed = christoffel(obj,ind,covIndArray)
            %christoffel i is the upper index and covIndArray is a 2-array
            %of indeces.
            %   Standard definition of the Christoffel symbols
            k = covIndArray(1);
            l = covIndArray(2);
            if ~obj.isComputedChristoffel
                obj.computeChristoffel();
            end

            christoffelIndexed = obj.grChristoffel(k,l,ind);
        end
        
        function riemannIndexed = riemann(obj, arg)
            % polymorphism in matlab is implemented via logic handling
            % with nargin and nargout
            if nargin == 4
                % assume R_{ijkl}
            end

            if nargin == 2
                % assume R^{i j}_{k l} or some other combo with non zero contravariant and covariant indeces
            end

            if nargin == 1
                % assume R^{i j k l}
            end
        end
                
    end

    methods (Access = private)
        function computeChristoffel(obj)
            %TODO: figure out a better way to compute this which is more efficient and uses symmetry
            obj.grChristoffel = sym('grStuff',[obj.grDimension,obj.grDimension,obj.grDimension]);
            for ind = 1:obj.grDimension
                for l = 1:obj.grDimension
                    for k = 1:obj.grDimension
                        chTerm = 0;
                        for m = 1:obj.grDimension
                            chTerm =chTerm + (1/2)* obj.grIMetric(ind,m)*(obj.grIMetricDerv(m,k,l) ...
                             + obj.grIMetricDerv(m,l,k) ...
                             - obj.grIMetricDerv(l,k,m) );
                        end                                                            
                        obj.grChristoffel(k,l,ind) = chTerm;
                    end
                end
            end
            obj.isComputedChristoffel = true;
        end


    end
end

