function padot = inverseDervProduct(invDervMetric,pa)
% pdot_a = p_b * p_c * dg^{bc}/dx^a 
%   Detailed explanation goes here
dim = length(pa);
padot = zeros(dim,1);
for a = 1 : dim
   padot(a) = -0.5 * dot(pa, invDervMetric(:,:,a) * pa); 
end
end

