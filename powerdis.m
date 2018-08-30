function [mainswitches,sideswitches]=powerdis(V0,Check)
I0=150;
disp("This is a 4*3 distribution network.");
%Maxcurrent=I0;
voltage=V0;
%disp("LOAD IS:")
load=randi(4,4,3) %#ok<NOPRT>
sideswitches=zeros(4,2);%Normally off two way switches.
%0 indicating the switch is off
%1 indicating the switch is on
%2 indicating the switch is normally connected to left and jump connected to  right and 3 for vice versa 
mainswitches=ones(4,3);%Normaaly on switches
n=size(mainswitches,1);
current=load./voltage*1000;
current1=current;
for i=1:n-1
current(i,1:end)=sum(mainswitches(i:end,1:end).*current(i:end,1:end));
end
disp("Intial current:-");
disp(current);
I=Check;
i=rem(I,4);
if i==0
    i=4;
end
j=fix(I/4)+1;
if I/4==fix(I/4)
    j=j-1;
end
fprintf('The probelm is at (%d,%d)\n',i,j);
mainswitches(i,j)=0;%indicating no current flows through it 
currentextra= (I0-current);%extra current a wire can carry when I0 is fixed
h=[];
for k=1:3
     if k~=j
         a=logical(current1(i:end,j)<=currentextra(1,k));
         h=[h,a];
    end
end
logic=h;
%disp(logic);
%disp(size(logic));
%not possilbe case
%==========================================================================
flag=0;
for k=1:3
    if k~=j
        if (~isempty(find(logic>0))) %#ok<EFIND>
            flag=1;
        end
    end
end
if flag==0
    disp(" power cannot be supplied to fault  area");
    for i=1:n-1
        current(i,1:end)=sum(mainswitches(i:end,1:end).*current1(i:end,1:end));
    end
    disp(current);
    return 
end
%==============================================================================
%Trivial case
for k=1:3
    if(k~=j)&&i~=1
        if currentextra(i,k)>current(i,j)
           if(i~=1)
               sideswitches(i-1,j-logical(j>k))=abs(k-j)+logical(abs(k-j)>1)*logical(j>k);%Because no loops should be in the diagram 
           end
               sideswitches(i,j-logical(j>k))=abs(k-j)+logical(logical(abs(k-j)>1)&logical(j>k));
            if(abs(k-j)>1)
               if i~=1 
                   sideswitches(i-1,k-logical(j<k))=abs(k-j)+logical(j<k);
               end
                sideswitches(i,k-logical(j<k))=abs(k-j)+logical(j<k);
            end
           current=Checkcurrent(mainswitches,sideswitches,i,j,current1,current,I0,currentextra);
           current=Checkcurrent1(current,current1,i,j,mainswitches,sideswitches,I0,currentextra);
           disp("Final Current in each bus is:-");
           disp(current);
           return;
        end
    end
end
%{
for k=1:3%Less efficient for small number of loads for a given feeder.High efficient for larger number of loads 
    for t=i+1:floor(n/2)
    if(k~=j)
        if currentextra(i,k)>current(t,j)
           if(i~=1)
               sideswitches(i-1,j-logical(j>k))=abs(k-j)+logical(abs(k-j)>1)*logical(j>k);%Because no loops should be in the diagram 
           end
               sideswitches(t,j-logical(j>k))=abs(k-j)+logical(logical(abs(k-j)>1)&logical(j>k));
            if(abs(k-j)>1)
               if i~=1 
                   sideswitches(i-1,k-logical(j<k))=abs(k-j)+logical(j<k);
               end
                sideswitches(t,k-logical(j<k))=abs(k-j)+logical(j<k);
            end
           current=Checkcurrent(mainswitches,sideswitches,i,j,current1,current,I0,currentextra);
           current=Checkcurrent1(current,current1,i,j,mainswitches,sideswitches,I0,currentextra);
           disp("Final Current in each bus is:-");
           disp(current);
           return;
        end
    end
    end
end
%}
%==========================================================================
%Non Trivial case     
disp("Non trivial case");
flag=0;
for k=1:3
    if(k~=j)
        if flag==0
            a=currentextra(1,k);
            flag=1;
        else
         b=currentextra(1,k);
        end
         if currentextra(1,k)>=current(i,j)
            sideswitches(i,j-logical(j>k))=1+logical(abs(k-j)>1)*logical(j>k);
            if(abs(k-j)>1)
                sideswitches(i,k-logical(j<k))=abs(k-j)+1;
            end
            current=Checkcurrent(mainswitches,sideswitches,i,j,current1,current,I0,currentextra);
            current=Checkcurrent1(current,current1,i,j,mainswitches,sideswitches,I0,currentextra);
            disp("Final Current in each bus is:-");
            disp(current);
            return;
         end
    end
end
%disp(a);
%disp(b);
%disp(logic);
logic3=supply(current1,logic,a,b,i,j);
disp("logic3:-") ;
disp(logic3);
[mainswitches,sideswitches]=supply1(logic3,i,j,mainswitches,sideswitches);
current=Checkcurrent(mainswitches,sideswitches,i,j,current1,current,I0,currentextra);
current=Checkcurrent1(current,current1,i,j,mainswitches,sideswitches,I0,currentextra);%Calculating final current for other buses
disp("Final Current in each bus is :-");
disp(current);
end
%============================================================================
%Checking by seeing currents in each wire 
%Calculating the last node current for each feeder
function current=Checkcurrent(mainswitches,sideswitches,i,j,current1,current,I0,currentextra)
for k=1:3
    current(end,k)=mainswitches(end,k)*current1(end,k);
    if(k<j)
        current(end,k)=current(end,k)+logical(sideswitches(end,k)==1)*current1(end,k+1);
        if abs(k-j)==2
           current(end,k)=current(end,k)+logical(sideswitches(end,k)==2)*current1(end,k+2);
        end
    elseif k>j
        current(end,k)=current(end,k)+logical(sideswitches(end,k-1)==1)*current1(end,k-1);        
        if abs(k-j)==2
           current(end,k)=current(end,k)+logical(sideswitches(end,k-1)==3)*current1(end,k-2);
        end  
    else 
        if j==1
           current(end,k)=current(end,k)+current1(end,k)*logical(sideswitches(end,k)~=0);
        elseif j==3
           current(end,k)=current(end,k)+current1(end,k)*logical(sideswitches(end,k-1)~=0);
        else
            current(end,k)=current(end,k)+current1(end,k)*(logical(sideswitches(end,k)~=0)||logical(sideswitches(end,k-1)~=0));
        end
    end
    current(end,k)=mainswitches(end,k)*current(end,k);
end
end

%==========================================================================
%{
function logic3=supply(logic)
logic3=zeros(size(logic,1),1);
x=size(logic,1);
logic2=or(logic(:,1),logic(:,2));
l=1;
    for t=1:x
        if(logic2(t)==1)
            p=find(logic(t,:)==1);%check here,Find gives linear indexing
            if sum(p==1)==1
            logic3(l)=1;
            else
                logic3(l)=2;
            end
        else
            logic3(l)=0;
        end
        l=l+1;
    end
    end

function logic3=supply(logic,current1,i,j,a,b)
    logic3=zeros(size(logic,1),1);
    x=size(logic,1);
    for c=1:(2^x)-1
        d=de2bi(c,x);
        combination(c)=sum(current1(i:end,j).*d(end-i+1));       
    end
    logic4=logical(combination<a+b);
    for c=1:2^x-1
        if(logic4(c)==1)
            d=de2bi(c,x);
           count(c)=sum(d);
        else
            count(c)=0;
        end
    end
   [~,I]= max(count);
   d=de2bi(I,x);
   l=1;
   hold=[];
    for c=1:x
        if d(c)~=1
            logic3(c)=0;
        end
        if d(c)==1
            if(current1(i+c-1,j)<=a&&current1(i+c-1,j)>b)
                logic3(c)=1;
                a=a-current1(i+c-1,j);
            elseif current1(i+c-1,j)<=b&&current1(i+c-1,j)>a
                logic3(c)=2;
                b=b-current1(i+c-1,j);
            elseif current1(i+c-1,j)<=b &&current1(i+c-1,j)<=a
                hold(l)=c;
                l=l+1;
            else
                logic3(c)=0;
            end
        end
    end
    y=size(hold,1);
    for t=1:y
        if(current1(hold(t)+i-1,j)<a)
             logic3(hold(t))=1;
                a=a-current1(i+hold(t)-1,j);
        elseif current1(hold(t)+i-1,j)<b
             logic3(hold(t))=2;
                b=b-current1(i+hold(t)-1,j);
        end
    end
end
%}
function logic3=supply(current1,logic,a,b,i,j)
x=size(logic,1);
logic3=zeros(x,1);
%current2=sort(current1);
for t=1:x
  
    if logic(t,1)==0&&logic(t,2)==0
        logic3(t,1:end)=0;
    elseif logic(t,1)==1 &&logic(t,2)==0
        logic3(t,1:end)=1;
    elseif logic(t,1)==0&&logic(t,2)==1
        logic3(t,1:end)=2;
    elseif logic(t,1)==1&&logic(t,2)==1
        logic3=[logic3,logic3];
        y=size(logic3,2);
        logic3(t,1:y/2)=1;
        logic3(t,(y/2)+1:end)=2;
    end
end
%disp("Medium logic3:-");
%disp(logic3);
logic4=logic3;
z=size(logic3,2);
for t=1:z
    c=a;
    d=b;
    for s=1:x
         if logic3(s,t)==1 && c>=current1(s+i-1,j)
           c=c-current1(s+i-1,j);
        elseif logic3(s,t)==2 && d>=current1(s+i-1,j)
            d=d-current1(s+i-1,j);
            logic3(s,t)=1;
        else 
            logic3(s,t)=0;
        end
    end
end
%disp(logic3);
for t=1:z
       b=sum(logic3,1);
end
%disp(b);
[~,I]=max(b);
%disp(I);
for t=1:x
    if logic3(t,I)==0
        logic4(t,I)=0;
    end
end
logic3=logic4(1:end,I);
end

%==========================================================================
function [mainswitches,sideswitches]= supply1(logic3,i,j,mainswitches,sideswitches)
   a=size(logic3,1);
    for t=1:a-1%Managing mainswitches such that they can get current from same sources or different
       
        if(logic3(t)==logic3(t+1)&&logic3(t)~=0)
            mainswitches(i+t,j)=1;
        else
            mainswitches(i+t,j)=0;
        end
    end
    for t=1:a
        if t>1
            if logic3(t)==logic3(t-1)
                continue;
            end
        end
        if(logic3(t)==1)
            if(j==1)
              sideswitches(i+t-1,j)=1;
            else 
                sideswitches(i+t-1,1)=abs(j-1);
                if(abs(j-1)>1)
                sideswitches(i+t-1,j-1)=3;  
                end
            end
        elseif (logic3(t)==2)
              if j==1
              sideswitches(i+t-1,j)=2;
              sideswitches(i+t-1,j+1)=3;
              else
                  sideswitches(i+t-1,j-logical(j>2))=1;
              end      
        end
    end
            
            
end
function current=Checkcurrent1(current,current1,i,j,mainswitches,sideswitches,I0,currentextra)
 n=size(mainswitches,1);
 for a=1:n-1
     l=n-a;
     for k=1:3
            if (mainswitches(l,k)==0)
                current(l,k)=0;
                continue;
            end
            current(l,k)=mainswitches(l,k)*(current1(l,k)+mainswitches(l+1,k)*current(l+1,k));
            if(k<j)
                current(l,k)=current(l,k)+logical(sideswitches(l,k)==1)*(mainswitches(l+1,k+1)*current(l+1,k+1)+current1(l,k+1))*logical(i<=l);
                if abs(k-j)==2
                    current(l,k)=current(l,k)+logical(sideswitches(l,k)==2)*(mainswitches(l+1,k+2)*current(l+1,k+2)+current1(l,k+2))*logical(i<=l);
                end
                if i>l&&sideswitches(l,k)~=0
                    current(l,k)=current(l,k)-current(l+1,k)+I0-currentextra(l+1,k);
                end
            elseif k>j
                current(l,k)=current(l,k)+logical(sideswitches(l,k-1)==1)*logical(i<=l)*(mainswitches(l+1,k-1)*current(l+1,k-1)+current1(l,k-1));        
                if abs(k-j)==2
                    current(l,k)=current(l,k)+logical(sideswitches(l,k-1)==3)*logical(i<=l)*(mainswitches(l+1,k-2)*current(l+1,k-2)+current1(l,k-2));
                end  
                 if i>l&&sideswitches(l,k-1)~=0
                    current(l,k)=current(l,k)-current(l+1,k)+I0-currentextra(l+1,k);
                end
            else
                 if j==1
                      current(l,k)=current(l,k)+(current(l+1,k+1)-I0+currentextra(l+1,k+1))*logical(sideswitches(l,k)==1)*logical(i>l);
                      if logical(sideswitches(l,k)==2)
                          current(l,k)=current(l,k)+(current(l+1,k+2)-I0+currentextra(l+1,k+2))*logical(sideswitches(l,k)==2)*logical(i>l);
                      end
                 elseif j==3
                      current(l,k)=current(l,k)+(current(l+1,k-1)-I0+currentextra(l+1,k-1))*logical(sideswitches(l,k-1)==1)*logical(i>l);
                      if logical(sideswitches(l,k-1)==3)
                          current(l,k)=current(l,k)+(current(l+1,k-2)-I0+currentextra(l+1,k-2))*logical(sideswitches(l,k-1)==3)*logical(i>l);
                      end
                 else
                      current(l,k)=current(l,k)+logical(sideswitches(l,k)==1)*logical(i>l)*(current(l+1,k+1)-I0+currentextra(l+1,k+1));
                      if logical(sideswitches(l,k-1)==1)
                          current(l,k)=current(l,k)+logical(i>l)*(current(l+1,k-1)-I0+currentextra(l+1,k-1));
                      end
                 end
            end
     end
 end
end
     



