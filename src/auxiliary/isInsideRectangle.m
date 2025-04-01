function d = isInsideRectangle(P,X1,X2)
d = (P(:,1)>=X1(1)) .* (P(:,2)>=X1(2)) .* (P(:,1)<=X2(1)) .* (P(:,2)<=X2(2));
