classdef gr_calc
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        grMetric;
        grCoordinates;
        grDimension;
        grIMetric;
        grMetricDerv = {};
        grIMetricDerv = {};
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
        
        function christoffelIndexed = christoffel(obj,i,covIndArray)
            %christoffel i is the upper index and covIndArray is a 2-array
            %of indeces.
            %   Standard definition of the Christoffel symbols
            k = covIndArray(1);
            l = covIndArray(2);
            gdervl = obj.grMetricDerv(:,:,l);
            gdervk = obj.grIMetricDerv(:,:,k);
            christoffelIndexed = 0;
            for m = 1:obj.grDimension
                christoffelIndexed =christoffelIndexed + (1/2)* obj.grIMetric(i,m)*(gdervl(m,k) + gdervk(m,l) - obj.grIMetricDerv(l,k,m) );
            end
        end
        
        function riemannIndexed = riemann(obj, arg)
            % polymorphism in matlab is implemented via logic handling
            % with nargin and nargout
            if nargin == 4
                % assume R_{ijkl}
            end
        end
                
    end
end

