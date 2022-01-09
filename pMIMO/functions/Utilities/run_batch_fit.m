%cd '/shared/persisted/';
%disp(pwd);
function F = run_batch_fit(node_fit,folder)
disp(node_fit)
F=node_fit.run(folder);