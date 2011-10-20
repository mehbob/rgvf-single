function [value,pos] = get_single_line(ori, ptStart, ptEnd)

assert(isa(ori,'double'),...
  'image input when computing the rad edge map is supposed to be double type');

y1 = ptStart(1);
x1 = ptStart(2);

y2 = ptEnd(1);
x2 = ptEnd(2);

d = sqrt(double((x1-x2)^2 + (y1-y2)^2));
d = round(d);

value = zeros(d+1, 1);
pos = zeros(d+1, 2);

value(1) = ori(y1,x1);
value(d+1) = ori(y2,x2);

pos(1,1)=y1;
pos(1,2)=x1;
pos(d+1,1)=y2;
pos(d+1,2)=x2;

for i = 2 : d
    xtmp = (x1*(d + 1 - i) + x2 *( i - 1))/d;
    ytmp = (y1*(d + 1 - i) + y2 *( i - 1))/d;
    
    pos(i,1) = ytmp;
    pos(i,2) = xtmp;
    
    value(i) = interp2(ori,xtmp,ytmp,'*linear');
end