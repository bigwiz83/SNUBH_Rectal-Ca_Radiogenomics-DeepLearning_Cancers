function data = matRead(filename)

inp = load(filename);
f = fields(inp);
data = inp.(f{1});
data = squeeze(data);
data = single(data);
Inputsize = [160 160 80];
data = imresize3(data, Inputsize, 'Method','linear');
