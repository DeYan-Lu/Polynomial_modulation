seq_length = 512; 
time_halfbandwidth = 2.5;
num_seq = 2*(3.5)-1;
[dps_seq,lambda] = dpss(seq_length,time_halfbandwidth,num_seq);

plot(dps_seq)
title('Slepian Sequences, N = 512, NW = 2.5')
axis([0 512 -0.15 0.15])
legend('1st','2nd','3rd','4th','5th','6th')
