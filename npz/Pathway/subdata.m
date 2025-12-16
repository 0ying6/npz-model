 clear,clc

fprintf("读取数据�?...");
x1 = load ('x1_trace.csv');
x3 = load ('x3_trace.csv');
x5 = load ('x5_trace.csv');
fprintf("已完成\n");
index=[330:336,690:695,1205:1210,1251:1255];
x1n=x1(index,:);
x3n=x3(index,:); 
x5n=x5(index,:);
% csvwrite('u_opt.csv',x1n)
% csvwrite('v_opt.csv',x1n)
save('subdatau_opt.mat','x1n')
save('subdatav_opt.mat','x3n')
save('subdataw_opt.mat','x5n')