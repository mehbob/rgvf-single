function varargout = maingui(varargin)
% MAINGUI MATLAB code for maingui.fig
%      MAINGUI, by itself, creates a new MAINGUI or raises the existing
%      singleton*.
%
%      H = MAINGUI returns the handle to a new MAINGUI or the handle to
%      the existing singleton*.
%
%      MAINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAINGUI.M with the given input arguments.
%
%      MAINGUI('Property','Value',...) creates a new MAINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before maingui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to maingui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help maingui

% Last Modified by GUIDE v2.5 15-Oct-2011 10:51:54

%global declaration
global img_rgb;
global img_gray;
global img_denoised;
global img_roughseg;
global mat_edgemap;
global mat_gvf;


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @maingui_OpeningFcn, ...
                   'gui_OutputFcn',  @maingui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before maingui is made visible.
function maingui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to maingui (see VARARGIN)

% Choose default command line output for maingui
handles.output = hObject;

pos=get(gcf,'Position');
set(gcf,'Position',[pos(1:2) 825 600]);
movegui(gcf, 'center');

%set the default values
handles.option_rgb2gray = 2;

%denoise
handles.denoise_method = 1;

handles.mean_neighbour_size = 3;
handles.median_neighbour_size = 3; % 3*3 window
handles.gaussian_hsize = 5;
handles.gaussian_sigma = 0.5;
handles.nlmeans_patch_size = 3;
handles.nlmeans_neighbour_size = 3;
handles.nlmeans_h = 20;
handles.bilateral_sigma_d = 3.0;
handles.bilateral_sigma_r = 25.0;
handles.bilateral_window_size = 3;
handles.ad_k = 0.05;
handles.ad_iter = 5;

%rough seg
handles.area_th_nu = 0.005;
handles.area_th_cyto = 0.2;
handles.rough_segmethod = 1; %Otsu


%rough seg
handles.option_comp_em_method = 2; %default: radiating
handles.centroid_sel_mode = 1; %1 for input the index and 2 for manually click 

handles.usersetcentroid = [0, 0];
handles.centroid_index = 1;
handles.stackrefine_theta = 2;
handles.positivesuppress_ratio = 0.1;
handles.ispositive_suppress = 1;
handles.isstack_refine = 1;

%gvf
handles.gvf_mu = 0.1;
handles.comp_gvf_iter = 100;
handles.is_gvf_normalize = 1;


%snake
handles.snake_dmax = 2;
handles.snake_dmin = 0.5;
handles.snake_alpha = 1;
handles.snake_beta = 1;
handles.snake_gamma = 1;
handles.snake_kappa = 1;
handles.snake_iter = 80;

handles.deform_nu_contour_mode = 1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes maingui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = maingui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% --- Executes during object creation, after setting all properties.
function axes_ori_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_ori (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');


% Hint: place code in OpeningFcn to populate axes_ori


% --- Executes during object creation, after setting all properties.
function axes_denoised_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_denoised (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');

% Hint: place code in OpeningFcn to populate axes_denoised


% --- Executes during object creation, after setting all properties.
function axes_em_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_em (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');

% Hint: place code in OpeningFcn to populate axes_em


% --- Executes during object creation, after setting all properties.
function axes_gvf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes_gvf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

set(hObject,'xTick',[]);
set(hObject,'ytick',[]);
set(hObject,'box','on');

% Hint: place code in OpeningFcn to populate axes_gvf


% --- Executes on button press in pb_back2grayscalepanel.
function pb_back2grayscalepanel_Callback(hObject, eventdata, handles)

% hObject    handle to pb_back2grayscalepanel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes_denoised);
cla;
set(handles.uipanel_denoising, 'Visible', 'Off');
set(handles.uipanel_rgb2grayscale, 'Visible', 'On');

% --- Executes on button press in pb_denoise.
function pb_denoise_Callback(hObject, eventdata, handles)
% hObject    handle to pb_denoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


global img_gray;
global img_denoised;

img_gray = BoundMirrorExpand(img_gray);
img_gray = img_gray / 255.0;

switch handles.denoise_method
    case 1
        %denoise using mean filter
        img_denoised = filter2(fspecial('average',...
                                        handles.mean_neighbour_size),...
                               img_gray);
        %since filter2 return type double, so normalize the img to [0,1]
        
    case 2
        %denoise using median filter
        img_denoised = medfilt2(img_gray, [handles.median_neighbour_size,...
                                           handles.median_neighbour_size]);
        
        
    case 3
        %denoise using gaussian filter
        stru_filt = fspecial('gaussian', handles.gaussian_hsize, ...
                        handles.gaussian_sigma);
        img_denoised = filter2(stru_filt, img_gray);
        
        
    case 4
        %denoise using anisotropic diffusion
        img_denoised = anisodiff(img_gray, handles.ad_iter, ...
                                    handles.ad_k, 0.25, 1);
         
    case 5
        %denoise using bilteral filter
         img_denoised  = bilateralfilter(img_gray,...
                                        [handles.bilateral_window_size,...
                                         handles.bilateral_window_size],...
                                         handles.bilateral_sigma_d,...
                                         handles.bilateral_sigma_r);
         
    case 6
        %denoise using nonlocal means filter
         Options.kernelratio = handles.nlmeans_patch_size;
         Options.windowratio = handles.nlmeans_neighbour_size;
         Options.filterstrength =  1.0 / handles.nlmeans_h;
         img_denoised = NLMF(img_gray,Options);
         
    otherwise
end

img_gray = BoundMirrorShrink(img_gray);

img_denoised = BoundMirrorShrink(img_denoised);

img_denoised = im2uint8(img_denoised);
img_gray = im2uint8(img_gray);

axes(handles.axes_denoised);
cla;
imshow(img_denoised);
set(handles.uipanel_denoising, 'Visible', 'Off');
set(handles.uipanel_roughsegmentation, 'Visible', 'On');
guidata(hObject, handles);





% --- Executes on selection change in popup_denoisingmethod.
function popup_denoisingmethod_Callback(hObject, eventdata, handles)
% hObject    handle to popup_denoisingmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_denoisingmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_denoisingmethod

handles = make_para_uipanels_invisible(handles);
switch get(handles.popup_denoisingmethod,'Value')
    case 1
        %make the corresponding selection panel visible
        handles.denoise_method = 1;
        set(handles.uipanel_mean_paras, 'Visible', 'On');
    case 2
        handles.denoise_method = 2;
        set(handles.uipanel_median_paras, 'Visible', 'On');
    case 3
        handles.denoise_method = 3;
        set(handles.uipanel_gaussian_paras, 'Visible', 'On');
    case 4
        handles.denoise_method = 4;
        set(handles.uipanel_ad_paras, 'Visible', 'On');
    case 5
        handles.denoise_method = 5;
        set(handles.uipanel_bilateral_paras, 'Visible', 'On');
    case 6
        handles.denoise_method = 6;
        set(handles.uipanel_nlmeans_paras, 'Visible', 'On');
    otherwise
end
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function popup_denoisingmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_denoisingmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_back2denoise.
function pb_back2denoise_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back2denoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_gray;
set(handles.uipanel_denoising, 'Visible', 'On');
set(handles.uipanel_roughsegmentation, 'Visible', 'Off');
axes(handles.axes_denoised);
imshow(img_gray);

% --- Executes on button press in pb_roughsegment.
function pb_roughsegment_Callback(hObject, eventdata, handles)
% hObject    handle to pb_roughsegment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_denoised;
global img_roughseg;
global centroids;
global boundaries_cyto;
global boundaries_nu;



switch handles.rough_segmethod
    case 1 % Otsu
        [img_roughseg,threshold_h, threshold_l] = otsu_1d3c(img_denoised);
    case 2 % Spatial kmeans
         img_roughseg = kmeans_spatial(img_denoised);
    otherwise
end


axes(handles.axes_denoised);
cla;
imshow(img_roughseg);

%get the initial contours and centroids
img_roughseg_reverse = 255 - img_roughseg;
thresholds = [0;99;199];
mask_cyto = getmask('cyto',img_roughseg_reverse,thresholds, img_denoised, ...
    handles.area_th_cyto); 
[mask_nu,props_nu] = getmask('nu',img_roughseg_reverse,...
                                  thresholds, img_denoised,...
                                  handles.area_th_nu, 0.7);

%select the largest cytoplam
tabu = tabulate(mask_cyto(:));
index1 = find(tabu(:,1)~=0); % tabu(:,1) == 0 means bkground
stat = tabu(index1,:);
stat = sortrows(stat,-2);

mask_cyto(mask_cyto ~= stat(1,1)) = 0;
mask_cyto = bwfill(mask_cyto, 'holes', 8);
boundaries_cyto = bwboundaries(mask_cyto);
num_cyto = size(boundaries_cyto,1);
assert(num_cyto == 1);


%deal with the nucleus
mask_nu(mask_cyto == 0) = 0;



tabu = tabulate(mask_nu(:));
index1 = find(tabu(:,1)~=0); % tabu(:,1) == 0 means bkground
stat = tabu(index1,:);
stat = sortrows(stat,-2);

num_nu = size(stat,1);
centroids = zeros(num_nu, 2);

for i = 1 : num_nu
    mat_help = mask_nu;
    mat_help(mat_help ~= stat(i,1)) = 0;
    
    centroids(i,1) = props_nu(stat(i,1),1).WeightedCentroid(1);
    centroids(i,2) = props_nu(stat(i,1),1).WeightedCentroid(2);
    
    tmp = bwboundaries(mat_help);
    boundaries_nu{i} = tmp{1,1};
end

%display the initial cyto and nu in the right down axes
axes(handles.axes_gvf);
cla;
imshow(img_denoised);
hold on;
plot(boundaries_cyto{1,1}(:,2),boundaries_cyto{1,1}(:,1),...
'Marker','*','LineWidth',1,'Color',[0 1 0]);

for i = 1 :  num_nu
    
    plot(boundaries_nu{i}(:,2),boundaries_nu{i}(:,1),...
        'Marker','*','LineWidth',1,'Color',[0 0 1]);
    plot(centroids(i,1),centroids(i,2),...
        'Marker','o','LineWidth',1,'Color',[1 0 0]);
    text(centroids(i,1),centroids(i,2),num2str(i),'LineWidth',1, 'Color',...
        [1,1,1],'FontSize',30);
 
end

set(handles.uipanel_roughsegmentation, 'Visible', 'Off');
set(handles.uipanel_compute_edgemap, 'Visible', 'On');

 set(handles.edit_centroid_index, 'String', ...
        num2str(handles.centroid_index));
guidata(hObject, handles);


% --- Executes on selection change in popupmenu_roughseg_method.
function popupmenu_roughseg_method_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_roughseg_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_roughseg_method contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_roughseg_method
handles.rough_segmethod = get(handles.popupmenu_roughseg_method,'Value');
 
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_roughseg_method_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_roughseg_method (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_back2roughseg.
function pb_back2roughseg_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back2roughseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_denoised;
set(handles.uipanel_roughsegmentation, 'Visible', 'On');
set(handles.uipanel_compute_edgemap, 'Visible', 'Off');
axes(handles.axes_gvf);
cla;
axes(handles.axes_denoised);
cla;
imshow(img_denoised);


% --- Executes on button press in pb_compute_edgemap.
function pb_compute_edgemap_Callback(hObject, eventdata, handles)
% hObject    handle to pb_compute_edgemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_edgemap;
global img_denoised;
global centroids;
global boundaries_cyto;
global boundaries_nu;

switch handles.option_comp_em_method
    case 1 % Sobel
        img_edgemap = sobelgradient(img_denoised);
    case 2
        %radiating edge map computation
        %get the position of centroid
        if(handles.centroid_sel_mode == 2)
            centroid = handles.usersetcentroid;
        else if (handles.centroid_sel_mode == 1)
                centroid = centroids(handles.centroid_index, :);
            end
        end
        Options.ispositive_suppress = handles.ispositive_suppress;
        Options.isstack_refine = handles.isstack_refine;
        Options.positivesuppress_ratio = handles.positivesuppress_ratio;
        Options.stackrefine_theta = handles.stackrefine_theta;
        img_edgemap = compute_rad_edgemap(double(img_denoised), centroid,...
                                        Options);
    otherwise
end
axes(handles.axes_em);
cla;
img_edgemap = im2uint8(img_edgemap / 255.0);
imshow(img_edgemap);
set(handles.uipanel_compute_edgemap, 'Visible', 'Off');
set(handles.uipanel_compute_gvf, 'Visible', 'On');
guidata(hObject, handles);

% --- Executes on button press in rb_sobel_em.
function rb_sobel_em_Callback(hObject, eventdata, handles)
% hObject    handle to rb_sobel_em (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_sobel_em


% --- Executes on button press in pb_setcentroid.
function pb_setcentroid_Callback(hObject, eventdata, handles)
% hObject    handle to pb_setcentroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% popup a big pic and let the user select the centroid
global img_gray;
h1 = figure;
imshow(img_gray);
handles.usersetcentroid = ginput(1);
close(h1);
guidata(hObject, handles);



% --- Executes on button press in cb_positivesuppress.
function cb_positivesuppress_Callback(hObject, eventdata, handles)
% hObject    handle to cb_positivesuppress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_positivesuppress

checkboxStatus = get(hObject,'Value');
if(checkboxStatus)
    %if box is checked, turn on positive suppress
    handles.ispositive_suppress = 1;
else
    %turn off
    handles.ispositive_suppress = 0;
end

guidata(hObject, handles);

% --- Executes on button press in cb_stackrefine.
function cb_stackrefine_Callback(hObject, eventdata, handles)
% hObject    handle to cb_stackrefine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_stackrefine
checkboxStatus = get(hObject,'Value');
if(checkboxStatus)
    %if box is checked, turn on positive suppress
    handles.isstack_refine = 1;
else
    %turn off
    handles.isstack_refine = 0;
end

guidata(hObject, handles);

% --- Executes on button press in pb_back2computeem.
function pb_back2computeem_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back2computeem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.uipanel_compute_edgemap, 'Visible', 'On');
set(handles.uipanel_compute_gvf, 'Visible', 'Off');
set(handles.edit_centroid_index, 'String', ...
        num2str(handles.centroid_index));
axes(handles.axes_em);
cla;

% --- Executes on button press in pb_compute_gvf.
function pb_compute_gvf_Callback(hObject, eventdata, handles)
% hObject    handle to pb_compute_gvf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_edgemap;
global img_gray;
global gvf_u;
global gvf_v;
global boundaries_cyto;
global boundaries_nu;
global centroids;

[gvf_u,gvf_v]= GVF(double(img_edgemap),handles.gvf_mu, handles.comp_gvf_iter); 

 if handles.is_gvf_normalize
    mag = sqrt(gvf_u.*gvf_u + gvf_v.*gvf_v);
    gvf_u = gvf_u./(mag+1e-10); 
    gvf_v = gvf_v./(mag+1e-10); 
 end
 
[m,n]=size(img_edgemap);

ySpace=(1:m/48 :m);
xSpace=(1:n/48 :n);

qx=interp2(gvf_u,xSpace, ySpace');
qy=interp2(gvf_v,xSpace, ySpace');
figure;
quiver(xSpace,ySpace,qx,qy); 

set(handles.uipanel_compute_gvf, 'Visible', 'Off');
set(handles.uipanel_snake_deform, 'Visible', 'On');

%set the edit_nucleus_index_deform status
if (size(centroids,1)  <= 1 ||...
    handles.centroid_sel_mode == 1 )
    set(handles.edit_nucleus_index_deform, 'Enable', 'Off');
    set(handles.edit_nucleus_index_deform, 'String', ...
        num2str(handles.centroid_index));
else
    set(handles.edit_nucleus_index_deform, 'Enable', 'On');
    set(handles.edit_nucleus_index_deform, 'String', ...
        num2str(handles.centroid_index));
end





function edit_mu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mu as text
%        str2double(get(hObject,'String')) returns contents of edit_mu as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.gvf_mu; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.gvf_mu = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_mu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_iter_comp_gvf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_iter_comp_gvf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_iter_comp_gvf as text
%        str2double(get(hObject,'String')) returns contents of edit_iter_comp_gvf as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.comp_gvf_iter; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.comp_gvf_iter = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_iter_comp_gvf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_iter_comp_gvf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_is_gvf_normalize.
function cb_is_gvf_normalize_Callback(hObject, eventdata, handles)
% hObject    handle to cb_is_gvf_normalize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_is_gvf_normalize
checkboxStatus = get(hObject,'Value');
if(checkboxStatus)
    %if box is checked, turn on positive suppress
    handles.is_gvf_normalize = 1;
else
    %turn off
    handles.is_gvf_normalize = 0;
end

guidata(hObject, handles);

% --- Executes on button press in pb_back2compgvf.
function pb_back2compgvf_Callback(hObject, eventdata, handles)
% hObject    handle to pb_back2compgvf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_rgb;
set(handles.uipanel_compute_gvf, 'Visible', 'On');
set(handles.uipanel_snake_deform, 'Visible', 'Off');
axes(handles.axes_ori);
cla;
imshow(img_rgb);


% --- Executes on button press in pb_deform_nu.
function pb_deform_nu_Callback(hObject, eventdata, handles)
% hObject    handle to pb_deform_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boundaries_nu;
global gvf_u;
global gvf_v;

%get the initial contour
switch handles.deform_nu_contour_mode
    case 1
        x = boundaries_nu{handles.centroid_index}(:,2);
        y = boundaries_nu{handles.centroid_index}(:,1);
    case 2
        %user draw the contour manually
        
    otherwise
end


[x,y] = snakeinterp(x,...
                    y,...
                    handles.snake_dmax,...
                    handles.snake_dmin);
                
dmax = handles.snake_dmax;
dmin = handles.snake_dmin;
alpha = handles.snake_alpha;
beta = handles.snake_beta;
gamma = handles.snake_gamma;
kappa = handles.snake_kappa;
iter = handles.snake_iter;

axes(handles.axes_ori);
hold on;
HD=line(x,y);
set(HD,'Color','Red');
for i=1:ceil(iter/5),
   if i<=floor(iter/5) 
      [x,y] = snakedeform(x,y,alpha,beta,gamma,kappa,gvf_u,gvf_v,5);
   else 
      [x,y] = snakedeform(x,y,alpha,beta,gamma,kappa,gvf_u,gvf_v,...
                        iter-floor(iter/5)*5);
   end;
   [x,y] = snakeinterp(x,y,dmax,dmin);
   
   XS=[x; x(1)];
   YS=[y; y(1)];
   
   delete(HD);
   HD=line(XS,YS);
   set(HD,'Color','Red');

   pause(0.1);
   
end



function edit_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_alpha as text
%        str2double(get(hObject,'String')) returns contents of edit_alpha as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.snake_alpha; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.snake_alpha = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_beta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_beta as text
%        str2double(get(hObject,'String')) returns contents of edit_beta as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.snake_beta; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.snake_beta = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_beta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_beta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_gamma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_gamma as text
%        str2double(get(hObject,'String')) returns contents of edit_gamma as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.snake_gamma; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.snake_gamma = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_gamma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_gamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_kappa_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kappa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kappa as text
%        str2double(get(hObject,'String')) returns contents of edit_kappa as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.snake_kappa; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.snake_kappa = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_kappa_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kappa (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dmin_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dmin as text
%        str2double(get(hObject,'String')) returns contents of edit_dmin as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.snake_dmin; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.snake_dmin = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_dmin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dmin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_dmax_Callback(hObject, eventdata, handles)
% hObject    handle to edit_dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_dmax as text
%        str2double(get(hObject,'String')) returns contents of edit_dmax as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.snake_dmax; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.snake_dmax = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_dmax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_dmax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_iter_deform_Callback(hObject, eventdata, handles)
% hObject    handle to edit_iter_deform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_iter_deform as text
%        str2double(get(hObject,'String')) returns contents of edit_iter_deform as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.snake_iter; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.snake_iter = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_iter_deform_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_iter_deform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pb_reset.
function pb_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pb_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close all;
clear;
clc;
maingui;

% --- Executes on button press in pb_draw_nu_contour.
function pb_draw_nu_contour_Callback(hObject, eventdata, handles)
% hObject    handle to pb_draw_nu_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_deform_cyto.
function pb_deform_cyto_Callback(hObject, eventdata, handles)
% hObject    handle to pb_deform_cyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global boundaries_cyto;
global gvf_u;
global gvf_v;

[x,y] = snakeinterp(boundaries_cyto{1}(:,2),...
                    boundaries_cyto{1}(:,1),...
                    handles.snake_dmax,...
                    handles.snake_dmin);
                
dmax = handles.snake_dmax;
dmin = handles.snake_dmin;
alpha = handles.snake_alpha;
beta = handles.snake_beta;
gamma = handles.snake_gamma;
kappa = handles.snake_kappa;
iter = handles.snake_iter;

axes(handles.axes_ori);
hold on;
HD=line(x,y);
set(HD,'Color','Red');
for i=1:ceil(iter/5),
   if i<=floor(iter/5) 
      [x,y] = snakedeform(x,y,alpha,beta,gamma,kappa,gvf_u,gvf_v,5);
   else 
      [x,y] = snakedeform(x,y,alpha,beta,gamma,kappa,gvf_u,gvf_v,...
                        iter-floor(iter/5)*5);
   end;
   [x,y] = snakeinterp(x,y,dmax,dmin);
   
   XS=[x; x(1)];
   YS=[y; y(1)];
   
   delete(HD);
   HD=line(XS,YS);
   set(HD,'Color','Red');

   pause(0.1);
   
end

% --- Executes on button press in pb_draw_cyto_contour.
function pb_draw_cyto_contour_Callback(hObject, eventdata, handles)
% hObject    handle to pb_draw_cyto_contour (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pb_openfile.
function pb_openfile_Callback(hObject, eventdata, handles)
% hObject    handle to pb_openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global img_rgb;
[file_name, file_path] = uigetfile({'*.bmp;*.jpg;*.png;*.pgm',...
                    'Image Files (*.bmp, *.jpg, *.png, *.pgm)'; ...
                    '*.*', 'All Files (*.*)'}, ...
                    'Pick a file');
img_rgb = imread([file_path file_name]);
[m,n,z] = size(img_rgb);
if m > n
    img_rgb = imresize(img_rgb, [200, round(200*n /m)]);
else
    img_rgb = imresize(img_rgb, [round(200*m/n), 200]);
end

axes(handles.axes_ori);
cla;
imshow(img_rgb);
set(handles.pb_rgb2grayscale,'Enable','On');
set(handles.pb_openfile,'Enable','Off');



% --- Executes on button press in pb_rgb2grayscale.
function pb_rgb2grayscale_Callback(hObject, eventdata, handles)
% hObject    handle to pb_rgb2grayscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global img_rgb;
global img_gray;

%get the option about how to transform rgb to grayscale
switch handles.option_rgb2gray
    case 1
        img_gray = rgb2gray(img_rgb);
    case 2
        struct_rgb2lab = makecform('srgb2lab');
        img_lab = applycform(img_rgb, struct_rgb2lab);
        ch_l = im2double(img_lab(:,:,1));
        min_l = min(ch_l(:));
        max_l = max(ch_l(:));
        if min_l ~= max_l
            norm_ch_l =  (ch_l - min_l) / (max_l - min_l);
        end
        img_gray = im2uint8(norm_ch_l);
    otherwise
end
axes(handles.axes_denoised);
cla;
imshow(img_gray);

set(handles.uipanel_denoising, 'Visible', 'On');
handles = make_para_uipanels_invisible(handles);

set(handles.uipanel_rgb2grayscale, 'Visible', 'Off');

switch get(handles.popup_denoisingmethod,'Value')
    case 1
        set(handles.uipanel_mean_paras, 'Visible', 'On');
    case 2
        set(handles.uipanel_median_paras, 'Visible', 'On');
    case 3
        set(handles.uipanel_gaussian_paras, 'Visible', 'On');
    case 4
        set(handles.uipanel_ad_paras, 'Visible', 'On');
    case 5
        set(handles.uipanel_bilateral_paras, 'Visible', 'On');
    case 6
        set(handles.uipanel_nlmeans_paras, 'Visible', 'On');
    otherwise
end

guidata(hObject, handles);



function edit_help_rgb2gray_Callback(hObject, eventdata, handles)
% hObject    handle to edit_help_rgb2gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_help_rgb2gray as text
%        str2double(get(hObject,'String')) returns contents of edit_help_rgb2gray as a double


% --- Executes during object creation, after setting all properties.
function edit_help_rgb2gray_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_help_rgb2gray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_options_rgb2gray.
function uipanel_options_rgb2gray_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_options_rgb2gray 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rb_traditional_rgb2gray'
        handles.option_rgb2gray = 1;

    case 'rb_cielab_rgb2gray'
        handles.option_rgb2gray = 2;
    otherwise
       % Code for when there is no match.
end
guidata(hObject, handles);



function edit_help_denoise_Callback(hObject, eventdata, handles)
% hObject    handle to edit_help_denoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_help_denoise as text
%        str2double(get(hObject,'String')) returns contents of edit_help_denoise as a double


% --- Executes during object creation, after setting all properties.
function edit_help_denoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_help_denoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_help_roughseg_Callback(hObject, eventdata, handles)
% hObject    handle to edit_help_roughseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_help_roughseg as text
%        str2double(get(hObject,'String')) returns contents of edit_help_roughseg as a double


% --- Executes during object creation, after setting all properties.
function edit_help_roughseg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_help_roughseg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_help_compute_edgemap_Callback(hObject, eventdata, handles)
% hObject    handle to edit_help_compute_edgemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_help_compute_edgemap as text
%        str2double(get(hObject,'String')) returns contents of edit_help_compute_edgemap as a double


% --- Executes during object creation, after setting all properties.
function edit_help_compute_edgemap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_help_compute_edgemap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_help_compute_gvf_Callback(hObject, eventdata, handles)
% hObject    handle to edit_help_compute_gvf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_help_compute_gvf as text
%        str2double(get(hObject,'String')) returns contents of edit_help_compute_gvf as a double


% --- Executes during object creation, after setting all properties.
function edit_help_compute_gvf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_help_compute_gvf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_help_deform_Callback(hObject, eventdata, handles)
% hObject    handle to edit_help_deform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_help_deform as text
%        str2double(get(hObject,'String')) returns contents of edit_help_deform as a double


% --- Executes during object creation, after setting all properties.
function edit_help_deform_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_help_deform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function handles = make_para_uipanels_invisible(handles)
            
set(handles.uipanel_mean_paras, 'Visible', 'Off');
set(handles.uipanel_median_paras, 'Visible', 'Off');
set(handles.uipanel_gaussian_paras, 'Visible', 'Off');
set(handles.uipanel_ad_paras, 'Visible', 'Off');
set(handles.uipanel_bilateral_paras, 'Visible', 'Off');
set(handles.uipanel_nlmeans_paras, 'Visible', 'Off');


% --- Executes during object creation, after setting all properties.
function uipanel_denoising_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_denoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit_mean_neighbourhood_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mean_neighbourhood_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mean_neighbourhood_size as text
%        str2double(get(hObject,'String')) returns contents of edit_mean_neighbourhood_size as a double

%get the paras

num = str2double(get(hObject,'String'));
if isnan(num) || round(num) ~= num || num <= 0
    num = handles.mean_neighbour_size; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a positive interger number', 'Error')
end
handles.mean_neighbourhood_size = num;
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function edit_mean_neighbourhood_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mean_neighbourhood_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_median_neighbourhood_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_median_neighbourhood_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_median_neighbourhood_size as text
%        str2double(get(hObject,'String')) returns contents of edit_median_neighbourhood_size as a double

%get the paras

num = str2double(get(hObject,'String'));
if isnan(num) || round(num) ~= num || num <= 0
    num = handles.median_neighbour_size; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a positive interger number', 'Error')
end

handles.median_neighbourhood_size = num;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit_median_neighbourhood_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_median_neighbourhood_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_denoise_gaussian_hsize_Callback(hObject, eventdata, handles)
% hObject    handle to edit_denoise_gaussian_hsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_denoise_gaussian_hsize as text
%        str2double(get(hObject,'String')) returns contents of edit_denoise_gaussian_hsize as a double

num = str2double(get(hObject,'String'));
if isnan(num) || round(num) ~= num || num <= 0
    num = handles.gaussian_hsize; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a positive interger number', 'Error')
end

handles.gaussian_hsize= num;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_denoise_gaussian_hsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_denoise_gaussian_hsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_denoise_gaussian_sigma_Callback(hObject, eventdata, handles)
% hObject    handle to edit_denoise_gaussian_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_denoise_gaussian_sigma as text
%        str2double(get(hObject,'String')) returns contents of edit_denoise_gaussian_sigma as a double

num = str2double(get(hObject,'String'));
if isnan(num) || num < 0
    num = handles.gaussian_sigma; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a non-negative number', 'Error')
end

handles.gaussian_sigma= num;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_denoise_gaussian_sigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_denoise_gaussian_sigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nlmeans_patch_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nlmeans_patch_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nlmeans_patch_size as text
%        str2double(get(hObject,'String')) returns contents of edit_nlmeans_patch_size as a double

num = str2double(get(hObject,'String'));
if isnan(num) || round(num) ~= num || num <= 0
    num = handles.nlmeans_patch_size; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a positive interger number', 'Error')
end

handles.nlmeans_patch_size = num;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_nlmeans_patch_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nlmeans_patch_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nlmeans_neighbour_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nlmeans_neighbour_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nlmeans_neighbour_size as text
%        str2double(get(hObject,'String')) returns contents of edit_nlmeans_neighbour_size as a double


num = str2double(get(hObject,'String'));
if isnan(num) || round(num) ~= num || num <= 0
    num = handles.nlmeans_neighbour_size; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a positive interger number', 'Error')
end

handles.nlmeans_neighbour_size = num;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_nlmeans_neighbour_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nlmeans_neighbour_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_nlmeans_h_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nlmeans_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nlmeans_h as text
%        str2double(get(hObject,'String')) returns contents of edit_nlmeans_h as a double

num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.nlmeans_h; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.nlmeans_h = num;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_nlmeans_h_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nlmeans_h (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ad_k_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ad_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ad_k as text
%        str2double(get(hObject,'String')) returns contents of edit_ad_k as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.ad_k; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.ad_k = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_ad_k_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ad_k (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_ad_iter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ad_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ad_iter as text
%        str2double(get(hObject,'String')) returns contents of edit_ad_iter as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.ad_iter; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.ad_iter = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_ad_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ad_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function axes_ori_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_ori (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_bilateral_window_size_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bilateral_window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bilateral_window_size as text
%        str2double(get(hObject,'String')) returns contents of edit_bilateral_window_size as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.bilateral_window_size; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.bilateral_window_size = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_bilateral_window_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bilateral_window_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bilateral_sigma_d_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bilateral_sigma_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bilateral_sigma_d as text
%        str2double(get(hObject,'String')) returns contents of edit_bilateral_sigma_d as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.bilateral_sigma_d; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.bilateral_sigma_d = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_bilateral_sigma_d_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bilateral_sigma_d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bilateral_sigma_r_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bilateral_sigma_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bilateral_sigma_r as text
%        str2double(get(hObject,'String')) returns contents of edit_bilateral_sigma_r as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.bilateral_sigma_r; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.bilateral_sigma_r = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_bilateral_sigma_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bilateral_sigma_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uipanel_inside_option_comp_em_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel_inside_option_comp_em (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when selected object is changed in uipanel_inside_option_comp_em.
function uipanel_inside_option_comp_em_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_inside_option_comp_em 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rb_sobel_em'
        handles.option_comp_em_method = 1;
        set(handles.edit_positivesuppress_ratio,'Enable','Off');
        set(handles.edit_stackrefine_theta,'Enable','Off');
        set(handles.cb_positivesuppress,'Enable','Off');
        set(handles.cb_stackrefine,'Enable','Off');
        set(handles.edit_centroid_index,'Enable','Off');
        set(handles.rb_input_nu_index,'Enable','Off');
        set(handles.rb_input_centroid_manually,'Enable','Off');

    case 'rb_rad_em'
        handles.option_comp_em_method = 2;
        set(handles.edit_positivesuppress_ratio,'Enable','On');
        set(handles.edit_stackrefine_theta,'Enable','On');
        set(handles.cb_positivesuppress,'Enable','On');
        set(handles.cb_stackrefine,'Enable','On');
        set(handles.edit_centroid_index,'Enable','On');
        set(handles.rb_input_nu_index,'Enable','On');
        set(handles.rb_input_centroid_manually,'Enable','On');
    otherwise
       % Code for when there is no match.
end
guidata(hObject, handles);



function edit_area_th_cyto_Callback(hObject, eventdata, handles)
% hObject    handle to edit_area_th_cyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_area_th_cyto as text
%        str2double(get(hObject,'String')) returns contents of edit_area_th_cyto as a double
num = str2double(get(hObject,'String'));
if isnan(num) || num < 0 || num >= 1
    num = handles.area_th_cyto; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number in (0,1)', 'Error')
end

handles.area_th_cyto = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_area_th_cyto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_area_th_cyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_area_th_nu_Callback(hObject, eventdata, handles)
% hObject    handle to edit_area_th_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_area_th_nu as text
%        str2double(get(hObject,'String')) returns contents of edit_area_th_nu as a double
num = str2double(get(hObject,'String'));
if isnan(num) || num < 0 || num >= 1
    num = handles.area_th_nu; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number in (0,1)', 'Error')
end

handles.area_th_nu = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_area_th_nu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_area_th_nu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_centroid_index_Callback(hObject, eventdata, handles)
% hObject    handle to edit_centroid_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_centroid_index as text
%        str2double(get(hObject,'String')) returns contents of edit_centroid_index as a double
global centroids;
num = str2double(get(hObject,'String'));
if isnan(num) || round(num) ~= num || num < 1 || num > size(centroids,1)
    num = handles.centroid_index; %use last value instead
    set(hObject,'String',num);
    errordlg('check input please', 'Error')
end

handles.centroid_index = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_centroid_index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_centroid_index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_positivesuppress_ratio_Callback(hObject, eventdata, handles)
% hObject    handle to edit_positivesuppress_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_positivesuppress_ratio as text
%        str2double(get(hObject,'String')) returns contents of edit_positivesuppress_ratio as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.positivesuppress_ratio; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.positivesuppress_ratio = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_positivesuppress_ratio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_positivesuppress_ratio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_stackrefine_theta_Callback(hObject, eventdata, handles)
% hObject    handle to edit_stackrefine_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_stackrefine_theta as text
%        str2double(get(hObject,'String')) returns contents of edit_stackrefine_theta as a double
num = str2double(get(hObject,'String'));
if isnan(num) 
    num = handles.stackrefine_theta; %use last value instead
    set(hObject,'String',num);
    errordlg('Input must be a number', 'Error')
end

handles.stackrefine_theta = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_stackrefine_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_stackrefine_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_sel_centroid.
function uipanel_sel_centroid_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_sel_centroid 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rb_input_nu_index'
        handles.centroid_sel_mode = 1;
        set(handles.pb_setcentroid,'Enable','Off');
        set(handles.edit_centroid_index,'Enable','On');

    case 'rb_input_centroid_manually'
        handles.centroid_sel_mode = 2;
        set(handles.edit_centroid_index,'Enable','Off');
        set(handles.pb_setcentroid,'Enable','On');
    otherwise
       % Code for when there is no match.
end
guidata(hObject, handles);



function edit_nucleus_index_deform_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nucleus_index_deform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nucleus_index_deform as text
%        str2double(get(hObject,'String')) returns contents of edit_nucleus_index_deform as a double
global centroids;
num = str2double(get(hObject,'String'));
if isnan(num) || round(num)~num || num < 1 || num > size(centroids,1)
    num = handles.centroid_index; %use last value instead
    set(hObject,'String',num);
    errordlg('check input please', 'Error')
end

handles.handles.centroid_index = num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_nucleus_index_deform_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nucleus_index_deform (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cb_deform_draw_cyto.
function cb_deform_draw_cyto_Callback(hObject, eventdata, handles)
% hObject    handle to cb_deform_draw_cyto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cb_deform_draw_cyto


% --- Executes when selected object is changed in uipanel_deform_nu_contour.
function uipanel_deform_nu_contour_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_deform_nu_contour 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag')   % Get Tag of selected object
    case 'rb_deform_input_nu_index'
        handles.deform_nu_contour_mode = 1; %use index

    case 'rb_deform_draw_nu'
        handles.deform_nu_contour_mode = 2; %draw contour manually
      
    otherwise
       % Code for when there is no match.
end
guidata(hObject, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
