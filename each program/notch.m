function [notch1,notch2,notch3]=notch(nf1,nf2,nf3,D2,D4,fs2,r1,r2,r3)
%ƒmƒbƒ`1  --------------------------

notch1 = [1 -2*r1*cos(nf1*pi/180) r1*r1];
if r1 == 0
    notch1 = 1;
else
    notch1 = notch1/sum(notch1);
end;
save notch1.cof notch1 -ascii;

%ƒmƒbƒ`2  ----------------------------

notch2 = [1 -2*r2*cos(nf2*pi/180) r2*r2];
if r2 == 0
    notch2 = 1;
else
    notch2 = notch2/sum(notch2);
end;
save notch2.cof notch2 -ascii;
%ƒmƒbƒ`3  ---------------------------

notch3 = [1 -2*r3*cos(nf3*pi/180) r3*r3];
if r3 == 0
    notch3 = 1;
else
    notch3 = notch3/sum(notch3);
end;
save notch3.cof notch3 -ascii;


