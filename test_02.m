fs = 1e4;    %サンプリング周波数10kHz
t = 0:1/fs:1;    

x =  sin(2*pi*60*t);
fo = 60;   %フィルタの次数

d = designfilt('hilbertfir', ...  %応答（フィルタの）タイプ
       'FilterOrder',fo, ...      %フィルタの次数
       'TransitionWidth',400, ... %遷移幅の指定→わからない
       'SampleRate',fs);          %サンプリング周波数

freqz(d,1024,fs)   %フィルタの周波数応答、出力引数なしだとフィルタの周波数応答がプロット
                   %(有理伝達関数のフィルタ係数？、n点分割、サンプリング周波数)

xh = filter(d,x);  %フィルタ処理の部分(フィルタ、入力データ)
                   %入力データxをフィルタdで処理
grd = fo/2;        %フィルタの群遅延
   
y2 = x(1:end-grd) + 1j*xh(grd+1:end);   %解析信号の生成
t2 = t(1:end-grd);                      %(1:end-grd)でベクトルにしてる？

plot(t2,real(y2),t2,imag(y2))  %複数線のプロットt2を独立変数(?)として実部と虚部を同一平面に記載
xlim([0.01 0.03])              %t=0.01~0.03までプロット
legend('real','imaginary')     %座標軸への凡例追加
title('FIR Filter')

instf=1/(2*pi)*diff(unwrap(angle(y2)));

plot(t2(2:end),instf(1:end))
ylim([0 10])
xlabel('Time (s)')
ylabel('Frequency (Hz)')


%瞬時新幅
