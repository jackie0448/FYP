function output = NDT(ref_img, input_img,iterations,varargin)

%ref img is the img that will be use as the final output //aka original
%input img is the color char that we will be using  // target
%check if the images fulfil the dimensional requirement
%{
ref_img= imread('ref_img.PNG');
input_img= imread('input_img.PNG');
imshow(ref_img);
imshow(input_img);
%}
if ((ndims(ref_img)~=3) && (ndims(input_img)~=3))
    error('images must be 3 dimensions');
end


numvarargs = length(varargin);
if numvarargs > 1
    error('pdf_transfer:TooManyInputs', ...
        'requires at most 1 optional input');
end

optargs = {1};
optargs(1:numvarargs) = varargin;
[relaxation] = optargs{:};


%this tells how many 'pages' are there in the 3rd Dimension
z_dim= size(ref_img,3);

%prof also mentioned at each iteration, take 1 channel
%note: channel refers to 1 page in the 3rd Dimension


%based on paper, change all the N dim into a 1 dim matrix in 1 step, in this case,
%basically taking each page of the 3 dim and convert to 1 dim- to process 1
%channel at a time
% a 3xM matrix where M can be any number
for(i = 1: z_dim)
    new_ref_img(i,:)= reshape(ref_img(:,:,i),1,size(ref_img,1)*size(ref_img,2));
    new_input_img(i,:)= reshape(input_img(:,:,i),1,size(input_img,1)*size(input_img,2));
end

%creating a rotation matrix uxing random projections where X(kxN) = Rk×dXd×N
% in this case, create a 1x3 rotation matrix??

%R{1}= [1 2];
R{1} = [1 2 0; 0 1 0; 0 0 1; 2/3 2/3 -1/3; 2/3 -1/3 2/3; -1/3 2/3 2/3];
for i= 1: iterations
   
     R{i}= R{1} * orth(randn(3,3)); %need to use rand to randomise it, orth to make it orth normal
  
   %rotate the coordinate system
   new_ref_img_rot= R{i} * new_ref_img;
   new_input_img_rot= R{i} * new_input_img;
   %to store ref img final output
    ref_img_final= ones(size(new_ref_img));
    
    preset_pixel= 300; %preset can be any value, min 250. this is just an approx to get a generic dataset
    num_of_proj= size(R{i},1); % based on the researh paper "automated color transfer", proj is like the axe
   
   %for each projections, get marginals, match and map the transformation
   for j= 1: num_of_proj
       %value to make inverse matrix possible
       eps = 1e-6;
       %get the unique marginals to map the 2 diff images tgt
       low_range= min([new_ref_img_rot(j,:) new_input_img_rot(j,:)]) - eps;
       max_range= max([new_ref_img_rot(j,:) new_input_img_rot(j,:)]) + eps;
       
       %to sorta defined an initialised dataset that meets the requirements
       u= ( (0:(preset_pixel-1))/ (preset_pixel-1) )* (max_range - low_range) + low_range;
       
       %no_bins= max_range-low_range; % maybe wrong
        
       %projections of samples based on the marginals, using histogram
       new_ref_img_proj= hist(new_ref_img_rot(j,:),u);
       new_input_img_proj= hist(new_input_img_rot(j,:),u);
       
        %use 0:no_bins bc we want to interpolate values for all the bins
       no_bins= max(size(new_ref_img_proj));
      
       %part 1 of the 1D pdf transfer, getting the cumulative graph
       ref_img_cum= cumsum(new_ref_img_proj + eps);
       ref_img_cum= ref_img_cum/ref_img_cum(numel(ref_img_cum));

       input_img_cum= cumsum(new_input_img_proj + eps);
       input_img_cum= input_img_cum/input_img_cum(numel(input_img_cum));
        
       %based on N-dim paper, returns the tN transformation from the cumsum
       %inversion
       
      %inverse_PY= inv (input_img_cum);
       
       look_up_table= interp1(input_img_cum,0:no_bins-1,ref_img_cum,'linear');
        min_benchmark= input_img_cum(1);
         max_benchmark= input_img_cum(end);
%this for loop is to check for possibilities of the probability in ref
%exceeding that in the input img

       for k= 1: no_bins
           if(ref_img_cum(k) < min_benchmark)
               look_up_table(k)= 0;
               
           end
           
           if (ref_img_cum(k) > max_benchmark)
               look_up_table(k)= no_bins-1;
               
           end
       
       end
       %}
       %{
       look_up_table(ref_img_cum<=input_img_cum(1))= 0;
       look_up_table(ref_img_cum>=input_img_cum(end))= no_bins-1;
       if sum(isnan(look_up_table))>0
           error('colour_transfer:pdf_transfer:NaN', ...
              'pdf_transfer has generated NaN values');
       end
       %}
       
       
       %apply mapping, aka tN(xN)
       ref_img_final(j,:) = (interp1(u, look_up_table, new_ref_img_rot(j,:))-1)/(preset_pixel-1)*(max_range-low_range) + low_range;
   end
   %whats this for
   new_ref_img = relaxation * (R{i} \ (ref_img_final - new_ref_img_rot)) + new_ref_img;
   
end
output= ref_img;
for i=1:z_dim
    output(:,:,i)= reshape(new_ref_img(i,:), size(output,1),size(output,2));
end


%TODO get 1Dimension PDF transfter functon, map back and apply mapping

