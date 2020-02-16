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
        
        function riemannIndexed = riemann(obj, arg)
        %   Computes the riemann tensor for a variety of index configurations
        %   see the documentation.
        %
            if ~isComputedRiemann
                obj.computeRiemann();
            end

            if nargin == 4
                % assume R_{ijkl}
                iv = arg{1};
                jv = arg{2};
                kv = arg{3};
                lv = arg{4};
                riemannIndexed = obj.grRiemann(iv,jv,kv,lv);
            end

            if nargin == 2
                % assume R^{i j}_{k l} or some other combo with non zero contravariant and covariant indeces
                contraIndArray = arg{1};
                covarIndArray = arg{2};
                if length(contraIndArray) == 1
                    iv = contraIndArray(1);
                    jv = covarIndArray(1);
                    kv = covarIndArray(2);
                    lv = covarIndArray(3);
                    riemannIndexed = obj.grRiemannUpper(jv,kv,lv,iv);
                else
                    %TODO TODO: FIXME: implement mixed index
                end
            end

            if nargin == 1
                % assume R^{i j k l}
            end
        end

        function ricciIndexed = ricci(obj, alph,bet)
            if ~isComputedRicci
                obj.computeRicci();
            end
            ricciIndexed = obj.grRicci(alph,bet);
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

        function computeRiemann(obj)
            % compute and store both R_{i,j,k,l} and R^i_jkl in order to avoid recomputation
            if ~obj.isComputedChristoffel
                obj.computeChristoffel();
            end
            dim = obj.grDimension;
            obj.grRiemann = sym('grStuff',dim,dim,dim,dim);
            obj.grRiemannUpper = sym('grStuff',dim,dim,dim,dim);
            % TODO: get rid of all these for loops and use the Riemann symmetry relations to improve storage
            for bet = 1:dim
                for gamm = 1:dim
                    for delt = 1:n
                        for alph = 1:dim
                            upperRiemannVal = diff(obj.grChristoffel(bet,delt,alph),obj.grCoordinates(gamm) - ...
                            diff(obj.grChristoffel(bet,gamm,alph),obj.grCoordinates(delt));

                            for mu = 1:n
                                upperRiemannVal = upperRiemannVal + ...
                                (obj.grChristoffel(mu,gamm,alph) * obj.grChristoffel(bet,delt,mu) - ...
                                obj.grChristoffel(mu,delt,alph) * obj.grChristoffel(bet,gamm,mu));
                            end
                            obj.grRiemannUpper(bet,gamm,delt,alph) = upperRiemannVal;
                        end
                        for cind = 1:dim
                            lowerRiemannVal = obj.grMetric(cind,1)*obj.grRiemannUpper(bet,gamm,delt,1);
                            for sumInd = 2:dim
                                lowerRiemannVal = lowerRiemannVal + obj.grMetric(cind,sumInd)*obj.grRiemannUpper(bet,gamm,delt,sumInd)
                            end
                            obj.grRiemann(bet,gamm,delt,cind) = lowerRiemannVal;                
                        end
                    end                    
                end                
            end            

            obj.isComputedRiemann = true;
        end

        function computeRicci(obj)
            if ~isComputedRiemann
                obj.computeRiemann();
            end
            dim = obj.grDimension;
            obj.grRicci = sym('grStuff',dim,dim);
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


    end
end

