function [NODE_D,FRACT] = generateFract(X0, X1, nDiv)

NODE_D = [linspace(X0(1),X1(1),nDiv+1)',linspace(X0(2),X1(2),nDiv+1)'];

FRACT = zeros(nDiv,2);
FRACT(1,:) = [1 2];
for i = 2:nDiv
    FRACT(i,:) = FRACT(i-1,:) + [1 1];
end

end
