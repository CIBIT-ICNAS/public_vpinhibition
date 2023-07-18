function [t,p] = slopeComparison(x1,y1,n1,x2,y2,n2)
    % based on https://www.graphpad.com/support/faq/how-does-prism-compare-the-slopes-of-linear-regression-lines-if-there-are-three-or-more-lines-how-can-i-account-for-multiple-comparisons/
    
    b1 = (x1*y1') / (x1*x1');
    b2 = (x2*y2') / (x2*x2');
    
    resSS1 = y1*y1';
    resSS2 = y2*y2';
    
    resDF1 = n1 - 2;
    resDF2 = n2 - 2;

    v = resDF1 + resDF2;
    
    s2 = (resSS1 + resSS2) / (resDF1 + resDF2);
    
    sb1b2 = sqrt( (s2/(x1*x1')) + (s2/(x2*x2')) );
    
    t = (b1 - b2) / sb1b2;
    
    critt = tinv(1-0.05,v);
    
    if t >= critt
        p = 2 * (1-tcdf(abs(t),v));
        fprintf('slopes are different: t = %0.4f, p=%0.4f\n',t,p);        
    end

end