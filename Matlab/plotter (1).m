figure(1)
imagesc(0:0.1:12.4,0:0.1:12.4,h)
title('Time of arrival histogram')
xlabel('Idler time of arrival (ns)')
ylabel('Signal time of arrival (ns)')

a = 35+1;
b = 25+1;
c = 5+1;
d = 11+1;

figure(2)
subplot(2,1,1)
plot([0:0.1:4],h(11:51,[c,d]),'LineWidth',2)
title('Idler histogram')
xlabel('Signal time of arrival given Idler (ns)')
legend(['Idler = ' num2str((c-1)/10) 'ns'],['Idler = ' num2str((d-1)/10) 'ns'])
ylabel('Counts')



subplot(2,1,2)
title('Signals histogram')
plot([0:0.1:4],h([a,b],1:41)','LineWidth',2)
xlabel('Idler time of arrival given Signal (ns)')
legend(['Signal = ' num2str((a-1)/10) 'ns'],['Signal = ' num2str((b-1)/10) 'ns'])
ylabel('Counts')




a = 35+1;
b = 25+1;
c = 5+1;
d = 13+1;



s1 = h(a,d);
s2 = h(a,c);
i1 = h(b,c);
i2 = h(b,d);

[s1,s2,i1,i2]

signals = [-ones(1,s1) +ones(1,s2) +ones(1,i1) -ones(1,i2)];
idlers = [ones(1,s1) +ones(1,s2) -ones(1,i1) -ones(1,i2)];


corrcoef(signals,idlers)


