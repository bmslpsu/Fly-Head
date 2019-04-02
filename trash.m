close all
figure(1)
hold on
plot(PAT.Freq{1,1}(:,1),PAT.Mag{1,1}(:,1))
plot(PAT.Peak.Freq{1,1}(:,1),PAT.Peak.Mag{1,1}(:,1))
hold off

figure(2)
plot(BODE.Freq{1,1},BODE.Mag.HeadPat{1,1})