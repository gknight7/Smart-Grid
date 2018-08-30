function [mainswitches,sideswitches]=powerdis2(V0,Check)
I0=150;
disp("This is a 4*3distribution network.");
%Maxcurrent=I0;
voltage=V0;
%disp("LOAD IS:")
load=randi(4,4,3)
sideswitches=zeros(4,2);%Normally off two way switches.
%0 indicating the switch is off
%1 indicating the switch is on
%2 indicating the switch is normally connected to left and jump connected to  right and 3 for vice versa 
mainswitches=ones(4,3);%Normaaly on switches
n=size(mainswitches,1);
%current=zeros(4,3);
current=load./voltage*1000;
current1=current;
for i=1:n-1
current(i,1:end)=sum(mainswitches(i:end,1:end).*current(i:end,1:end));
end
%disp("Initial current is");
disp("Intial current:-");
disp(current);
%Check(2,1)=1;
%I=find(Check>0);
I=Check;
i=rem(I,4);
if i==0
    i=4;
end
j=fix(I/4)+1;
if I/4==fix(I/4)
    j=j-1;
end
mainswitches(i,j)=0;%indicating no current flows through it 
currentextra= (I0-current);%extra current a wire can carry when I0 is fixed
h=[];
for k=1:3
     if k~=j
         a=logical(current1(i:end,j)<currentextra(1,k));
         h=[h,a];
    end
end
logic=h;
size(logic);
%not possilbe case
%==========================================================================
flag=0;
for k=1:3
    if k~=j
        if (~isempty(find(logic>0)))
            flag=1;
        end
    end
end
if flag==0
    disp(" power cannot be supplied to fault  area");
    return 
end
%==============================================================================
%Trivial case
for k=1:3
    if(k~=j)
        if currentextra(i,k)>current(i,j)
            sideswitches(i-1,j-logical(j>k))=abs(k-j)+logical(abs(k-j)>1);%Because no loops should be in the diagram 
            sideswitches(i,j-logical(j>k))=abs(k-j)+logical(abs(k-j)>1);
            if(abs(k-j)>1)
                sideswitches(i-1,k-logical(j<k))=abs(k-j)+1;
                sideswitches(i,k-logical(j<k))=abs(k-j)+1;
            end
            return;
        end
    end
end
%==========================================================================
%Non Trivial case       
for k=1:3
    if(k~=j)
         if currentextra(1,k)>current(i,j)
            sideswitches(i,j-logical(j>k))=abs(k-j)+logical(abs(k-j)>1);
            if(abs(k-j)>1)
                sideswitches(i,k-logical(j<k))=abs(k-j)+1;
            end
            return;
         else 
           logic3=supply(logic); 
           [mainswitches,sideswitches]=supply1(logic3,i,j,mainswitches,sideswitches);
           
         end
    end
end
%============================================================================
%{
Checking by seeing currents in each wire 
a=size(mainswitches,1);
for i=1:3
    for t = 0:a-i
        
    effective_switch(i,1:end)=mainswitches(a-t,1:end)*(1||effective_switch(i,1:end));%A*(1||B*(1||C))
    end
end
for i=1:3
    current(i,1:end)=sum(effective_switch(i:end,1:end).*current1(i:end,1:end));
end

%}
disp('fuc');
end

%==========================================================================
function logic3=supply(logic)
logic3=zeros(size(logic,1),1);
x=size(logic,1);
logic2=or(logic(:,1),logic(:,2));
l=1;
    for t=1:x
        if(logic2(t)==1)
            p=find(logic(t,:)==1);%check here
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
%==========================================================================
function [mainswitches,sideswitches]= supply1(logic3,i,j,mainswitches,sideswitches)
   a=size(mainswitches,1);
    for t=i:a-i
        if(logic3(t)==logic3(t+1))
            mainswitches(t+1)=1;
        else
            mainswitches(t+1)=0;
        end
    end
    for t=i:a-i
        if(logic3(t)==1)
            if(j==1)
              sideswitches(i,j)=1;
            else 
                sideswitches(i,1)=abs(j-1);
                if(abs(j-1)>1)
                sideswitches(i,j-1)=3;  
                end
            end
        elseif (logic3(t)==2)
              if j==1
              sideswitches(i,j)=2;
              sideswitches(i,j+1)=3;
              else
                  sideswitches(i,j-logical(j>2))=1;
              end
        end
    end
            
            
end