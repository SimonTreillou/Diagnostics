h=linspace(1,10,20);
beta=h./180;
L=2*pi ./w_to_k(2*pi/6,h);
cond= (1+6.4 .*beta) .* h./L;

plot(h,cond)
hold on
plot(h,beta)


NP_XI=4
Lm=(372+NP_XI-1)/NP_XI;
padd_X=(Lm+2)/2-(Lm+1)/2;
disp(Lm+padd_X);

260
605
