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

cgif = menu("Done.", "Plot real-time", "Create GIF file");

if cgif == 1
	if isempty(findall(0,'Type','Figure'))
		figure()
	end
	
	clf;
	disp("Click on figure to plot..."); fflush(stdout);
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

	hold on

	for i = 1:n
		plot(data(:,2*i),data(:,2*i+1),"color", mycmap(i,:), "linewidth", 1);
	end

	hold off
end

if cgif == 2
	disp("Printing frames..."); fflush(stdout);
	set(0, 'defaultfigurevisible', 'off');
	mkdir output;
	for t = 1:size(normaldata,1)
		scatter(normaldata(t,[2:2:end]),normaldata(t,[3:2:end]),5,mycmap(1:n,:),"filled");
		grid on
		axis(baseaxis,"equal");
		title(["t = " num2str(normaldata(t,1),"%2.2f")]);
		filename=sprintf('./output/frame%d.png',t);
		print("-r0", filename);
	end
	set(0, 'defaultfigurevisible', 'on');
	disp("Done. Converting and creating gif..."); fflush(stdout);
	system("for i in $(ls ./output/*.png); do convert $i $i.gif; rm $i; done");
	system("gifsicle --delay=3 $(ls -v1 ./output/*.gif) > output.gif");
	system("rm -r ./output");
end

disp("Done!")