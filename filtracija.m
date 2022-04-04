function ctrlf = filtracija(i,ctrl)
CTRL(i)=ctrl;
if i==1
     ctrlf=0.4*CTRL(i);
elseif i==2
     ctrlf=0.4*CTRL(i)+0.3*CTRL(i-1);
elseif i==3
     ctrlf=0.4*CTRL(i)+0.3*CTRL(i-1)+0.2*CTRL(i-2);
else
     ctrlf=0.3*CTRL(i)+0.35*CTRL(i-1)+0.2*CTRL(i-2)+0.15*CTRL(i-3);
end
end

