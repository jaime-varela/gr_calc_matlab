classdef gr_calc  < handle
    %gr_calc a class to compute standard GR objects
    %   constructor: gr_calc(metric,coordinates)
    %   various methods to compute standard general relativity objects.
    

    % Todo: valid index function which throws if index values and types are invalid

    properties
        grMetric;
        grCoordinates;
        grDimension;
        grIMetric;
        grMetricDerv = {};
        grIMetricDerv = {};

        % usage \Gamma^i_{kl} = this.grChristoffel(k,l,i)
        grChristoffel = {};
        % usage R_{ijkl} = this.grRiemann(i,j,k,l)
        grRiemann = {};
        % usage R^i_{jkl} = this.grRiemannUpper(j,k,l,i)
        grRiemannUpper = {};
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
        
        function riemannIndexed = riemann(obj, varargin)
        %   Computes the riemann tensor for a variety of index configurations
        %   see the documentation.
        %
            if ~obj.isComputedRiemann
                obj.computeRiemann();
            end
            numArgs = nargin -1;
            if numArgs == 4
                % assume R_{ijkl}
                iv = varargin{1};
                jv = varargin{2};
                kv = varargin{3};
                lv = varargin{4};
                riemannIndexed = obj.grRiemann(iv,jv,kv,lv);
            end

            if numArgs == 2
                % assume R^{i}_{j k l} or some other combo with non zero contravariant and covariant indeces
                    iv = varargin{1};
                    covarIndArray = varargin{2};
                    jv = covarIndArray(1);
                    kv = covarIndArray(2);
                    lv = covarIndArray(3);
                    riemannIndexed = obj.grRiemannUpper(jv,kv,lv,iv);
            end

            if numArgs == 1
                % assume R^{i j k l} with an index array as input
                %TODO: implement the sum
            end
        end

        function ricciIndexed = ricci(obj, alph,bet)
            if ~obj.isComputedRicci
                obj.computeRicci();
            end
            ricciIndexed = obj.grRicci(alph,bet);
        end

        function rScalar = ricciScalar(obj)
            if ~obj.isComputedRScalar
                obj.computeRScalar();
            end
            rScalar = obj.grRscalar;            
        end

        function einsteinIndexed = einstein(obj,alph,bet)
            if ~obj.isComputedEinstein
                obj.computeEinstein();
            end

            einsteinIndexed = obj.grEinstein(alph,bet);            
        end

        function computeAllObjects(obj)
            % the following call should calculate and store all objects
            obj.einstein(0,0);
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
                            chTerm =chTerm + (1/2)* obj.grIMetric(ind,m)*(obj.grMetricDerv(k,m,l) ...
                             + obj.grMetricDerv(m,l,k) ...
                             - obj.grMetricDerv(l,k,m) );
                        end                                                            
                        obj.grChristoffel(k,l,ind) = chTerm;
                    end
                end
            end
            obj.isComputedChristoffel = true;
        end

        function computeRiemann(obj)
            % compute and store both R_{i,j,k,l} and R^i_jkl in order to avoid recomputation
            if ~obj.isComputedChristoffel
                obj.computeChristoffel();
            end
            dim = obj.grDimension;
            obj.grRiemann = sym('grStuff',[dim,dim,dim,dim]);
            obj.grRiemannUpper = sym('grStuff',[dim,dim,dim,dim]);
            % TODO: get rid of all these for loops and use the Riemann symmetry relations to improve storage
            for bet = 1:dim
                for gamm = 1:dim
                    for delt = 1:dim
                        for alph = 1:dim
                            upperRiemannVal = diff(obj.grChristoffel(bet,delt,alph),obj.grCoordinates(gamm)) - ...
                            diff(obj.grChristoffel(bet,gamm,alph),obj.grCoordinates(delt));

                            for mu = 1:dim
                                upperRiemannVal = upperRiemannVal + ...
                                (obj.grChristoffel(mu,gamm,alph) * obj.grChristoffel(bet,delt,mu) - ...
                                obj.grChristoffel(mu,delt,alph) * obj.grChristoffel(bet,gamm,mu));
                            end
                            obj.grRiemannUpper(bet,gamm,delt,alph) = upperRiemannVal;
                        end
                        for cind = 1:dim
                            lowerRiemannVal = obj.grMetric(cind,1)*obj.grRiemannUpper(bet,gamm,delt,1);
                            for sumInd = 2:dim
                                lowerRiemannVal = lowerRiemannVal + obj.grMetric(cind,sumInd)*obj.grRiemannUpper(bet,gamm,delt,sumInd);
                            end
                            obj.grRiemann(cind,bet,gamm,delt) = lowerRiemannVal;                
                        end
                    end                    
                end                
            end            

            obj.isComputedRiemann = true;
        end

        function computeRicci(obj)
            if ~obj.isComputedRiemann
                obj.computeRiemann();
            end
            dim = obj.grDimension;
            obj.grRicci = sym('grStuff',[dim,dim]);
            for alph = 1:dim
                for bet = 1:dim
                    ricciVal = obj.grRiemannUpper(alph,1,bet,1);
                    for sumInd = 2:dim
                        ricciVal = ricciVal + obj.grRiemannUpper(alph,sumInd,bet,sumInd);
                    end
                    obj.grRicci(alph,bet) = ricciVal;
                end
            end
            obj.isComputedRicci = true;            
        end

        function computeRScalar(obj)
            if ~obj.isComputedRicci
                obj.computeRicci();
            end
            obj.grRscalar = trace(obj.grIMetric * obj.grRicci);
            obj.isComputedRScalar = true;                        
        end

        function computeEinstein(obj)
            if ~obj.isComputedRScalar
                obj.computeRScalar();
            end

            obj.grEinstein = obj.grRicci - (1/2)*obj.grMetric * obj.grRscalar;            
        end


    end
end

