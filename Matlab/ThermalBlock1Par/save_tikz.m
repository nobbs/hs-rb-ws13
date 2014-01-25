form = '%.6g';

figure(1)
matlab2tikz('test1.tikz', 'encoding', 'utf8', 'height', '\figureheight', 'width', '\figurewidth', 'floatFormat', form)

figure(2)
matlab2tikz('test2.tikz', 'encoding', 'utf8', 'height', '\figureheight', 'width', '\figurewidth', 'floatFormat', form)

figure(3)
matlab2tikz('test3.tikz', 'encoding', 'utf8', 'height', '\figureheight', 'width', '\figurewidth', 'floatFormat', form)

figure(4)
matlab2tikz('test4.tikz', 'encoding', 'utf8', 'height', '\figureheight', 'width', '\figurewidth', 'floatFormat', form)
