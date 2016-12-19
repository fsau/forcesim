disp("Calculating..."); fflush(stdout);

fflag = system("gcc sim.c -O3 -lm -o sim && ./sim > orbit1.dat");
if fflag != 0
	return;
end

disp("Done. Processing..."); fflush(stdout);

fps = 30;

fh = fopen("./orbit1.dat");
raw = fscanf(fh, "%f");
fclose(fh);

n = raw(1); % Number of particles
raw = raw(2:end);

data = reshape(raw, [2*n+1, length(raw)/(2*n+1)])';

tmax = max(data(:,1));
normaldata = [(0:1/fps:tmax)', spline(data(:,1), data(:,2:end),0:1/fps:tmax)'];
datalength = size(normaldata,1);

baseaxis = [min(min(data(:,[2:2:end]))),max(max(data(:,[2:2:end]))),...
min(min(data(:,[3:2:end]))),max(max(data(:,[3:2:end])))];
mycmap = colormap(lines);

clf;
disp("Done. Click figure to plot..."); fflush(stdout);
waitforbuttonpress();

i = 1;
tic;
while i < datalength
	scatter(normaldata(i,[2:2:end]),normaldata(i,[3:2:end]),5,mycmap(1:n,:),"filled");
	grid on
	axis(baseaxis,"equal");
	title(["t = " num2str(normaldata(i,1),"%2.2f")]);
	pause(0.001);
	i = round(fps*toc);
end

disp("Done!")