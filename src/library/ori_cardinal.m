function matout = ori_cardinal(matin)

if size(matin,1)~=24
    matout = nan;
else
    p1 = [matin(19:24,:); matin(1:7,:)]; 
    p2 = matin(7:19,:); 
    matout = nan(size(p1)); 
    for ii = 1:size(matin,2)
        matout(:,ii) = circ_m([p1(:,ii) p2(:,ii)]'); 
    end
end