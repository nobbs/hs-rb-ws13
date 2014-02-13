% Path erweitern
addpath('./geometries');
addpath('./tools');

clear;

disp('Lade Geometrie des Definitionsgebiets.');
zweixeins;

disp('Bereite Offline-Phase vor.');
prepare_offline_stage;

disp('Beginne Offline-Phase.');
offline_stage;
