%------Task 1-------
function y_pred = lagrange_interpolation(x, x_points, y_points)
    y = 0;
    n = length(x_points);
    for i = 1:n
        p = 1;
        for j = 1:n
            if i ~= j
                p = p * (x - x_points(j))/(x_points(i) - x_points(j)); %Product notation formula
            end
        end
        y = y + p * y_points(i);
    end
    y_pred = y; 
end

function y_pred = cubic_spline (x, x_points, y_points)
    a_values = y_points;
    n = length(x_points) - 1; % number of subintervals 
    %Step 1
    h_values = zeros(1, n);
    for i = 1:n
        h_values(i) = x_points(i+1) - x_points(i);
    end
    %Step 2
    alpha_values = zeros(1, n+1);
    for i = 2:n
        alpha_values(i) = 3 * (a_values(i+1) - a_values(i))/h_values(i) - 3 * (a_values(i) - a_values(i-1))/h_values(i-1);
    end
    %Step 3 
    l_values = ones(1, n+1);
    m_values = zeros(1, n);
    z_values = zeros(1, n+1);
    %Step 4 
    for i = 2:n
        l_values(i) = 2 * (x_points(i+1) - x_points(i-1)) - h_values(i-1) * m_values(i-1);
        m_values(i) = h_values(i) / l_values(i);
        z_values(i) = (alpha_values(i) - h_values(i-1)*z_values(i-1))/(l_values(i));
    end
    %Step 5
    c_values = zeros(1, n+1);
    %Step 6 
    b_values = zeros(1, n);
    d_values = zeros(1, n);
    for j = n:-1:1
        c_values(j) = z_values(j) - m_values(j) * c_values(j+1);
        b_values(j) = (a_values(j+1) - a_values(j))/(h_values(j)) - h_values(j) * (c_values(j+1) + 2 * c_values(j))/3;
        d_values(j) = (c_values(j+1) - c_values(j))/(3 * h_values(j));
    end

     n = length(x_points); %Now set n to be the number of the points
    index = 0;
    for i = 2:n
        if x >= x_points(i-1) && x <= x_points(i)
            index = i-1;
            break;
        elseif x > x_points(n)
            index = n-1;
            break;
        elseif x < x_points(1)
            index = 1;
            break;
        end
    end

    y_pred = a_values(index) + b_values(index) * (x - x_points(index)) + c_values(index) * (x - x_points(index))^2 + d_values(index) * (x - x_points(index))^3;
end

% Inputs
year_values = [2019, 2020, 2021, 2022, 2023, 2024];
all_pop = [18395567, 18631779, 18879552, 19503159, 19776807, 20033842];
urban_pop = [10698208, 10938652, 11151376, 11991238, 12209896, 12451192];
rural_pop = [7697359, 7693127, 7728176, 7511921, 7556911, 7582650];

% Task 1
%Lagrange interpolation
disp("Lagrange interpolations results: ");
lpop_all = lagrange_interpolation(2024, year_values(2:5), all_pop(2:5));
lpop_all_abs = abs(all_pop(6) - lpop_all);
lpop_all_rel = round((lpop_all_abs/all_pop(6)), 5);

lpop_urban = lagrange_interpolation(2024, year_values(2:5), urban_pop(2:5));
lpop_urban_abs = abs(urban_pop(6) - lpop_urban);
lpop_urban_rel = round((lpop_urban_abs/urban_pop(6)), 5);

lpop_rural = lagrange_interpolation(2024, year_values(2:5), rural_pop(2:5));
lpop_rural_abs = abs(rural_pop(6) - lpop_rural);
lpop_rural_rel = round((lpop_rural_abs/rural_pop(6)), 5);

formattedStr = sprintf("Predicted population for year %d is %d (Abs. error: %i, Rel. error: %.3d)", 2024, lpop_all, lpop_all_abs, lpop_all_rel);
disp(formattedStr);
formattedStr = sprintf("Predicted urban population for year %d is %d (Abs. error: %i, Rel. error: %.3d)", 2024, lpop_urban, lpop_urban_abs, lpop_urban_rel);
disp(formattedStr);
formattedStr = sprintf("Predicted rural population for year %d is %d (Abs. error: %i, Rel. error: %.3d)", 2024, lpop_rural, lpop_rural_abs, lpop_rural_rel);
disp(formattedStr);


%For cubic spline interpolation
disp("Cubic spline interpolation results:");
spop_all = cubic_spline(2024, year_values(2:5), all_pop(2:5));
spop_all_abs = abs(all_pop(6) - spop_all);
spop_all_rel = round((spop_all_abs/all_pop(6)), 5);
spop_urban = cubic_spline(2024, year_values(2:5), urban_pop(2:5));
spop_urban_abs = abs(urban_pop(6) - spop_urban);
spop_urban_rel = round((spop_urban_abs/urban_pop(6)), 5);
spop_rural = cubic_spline(2024, year_values(2:5), rural_pop(2:5));
spop_rural_abs = abs(rural_pop(6) - spop_rural);
spop_rural_rel = round((spop_rural_abs/rural_pop(6)), 5);

formattedStr = sprintf("Predicted population for year %d is %d (Abs. error: %i, Rel. error: %.3d)", 2024, spop_all, spop_all_abs, spop_all_rel);
disp(formattedStr);
formattedStr = sprintf("Predicted urban population for year %d is %d (Abs. error: %i, Rel. error: %.3d)", 2024, spop_urban, spop_urban_abs, spop_urban_rel);
disp(formattedStr);
formattedStr = sprintf("Predicted rural population for year %d is %d (Abs. error: %i, Rel. error: %.3d)", 2024, spop_rural, spop_rural_abs, spop_rural_rel);
disp(formattedStr);

%Plotting for lagrange interpolation
all_modified = all_pop;
all_modified(6) = lpop_all;
urban_modified = urban_pop;
urban_modified(6) = lpop_urban;
rural_modified = rural_pop;
rural_modified(6) = lpop_rural;

figure
plot(year_values, all_pop, '-', 'LineWidth', 2)  
hold on                           
plot(year_values, all_modified, ':', 'LineWidth', 2)
plot(year_values, urban_pop, '-', 'LineWidth', 2)  
plot(year_values, urban_modified, ':', 'LineWidth', 2)
plot(year_values, rural_pop, '-', 'LineWidth', 2)  
plot(year_values, rural_modified, ':', 'LineWidth', 2)

xlabel('Year')
ylabel('Population')
title('Lagrange interpolation')
legend('Actual data (all)', 'Lagrange interpolation (all)', 'Actual data (urban)', 'Lagrange interpolation (urban)', 'Actual data (rural)', 'Lagrange interpolation (rural)')
grid on
hold off

%Plotting for cubic spline
all_modified = all_pop;
all_modified(6) = spop_all;
urban_modified = urban_pop;
urban_modified(6) = spop_urban;
rural_modified = rural_pop;
rural_modified(6) = spop_rural;

figure
plot(year_values, all_pop, '-', 'LineWidth', 2)  
hold on                           
plot(year_values, all_modified, ':', 'LineWidth', 2)
plot(year_values, urban_pop, '-', 'LineWidth', 2)  
plot(year_values, urban_modified, ':', 'LineWidth', 2)
plot(year_values, rural_pop, '-', 'LineWidth', 2)  
plot(year_values, rural_modified, ':', 'LineWidth', 2)

xlabel('Year')
ylabel('Population')
title('Cubic spline')
legend('Actual data (all)', 'Cubic spline (all)', 'Actual data (urban)', 'Cubic spline (urban)', 'Actual data (rural)', 'Cubic spline (rural)')
grid on
hold off

disp(" ")

% --------- Task 4 ----------
M = 20200000;
acc = 10;
x = 2019;
h = 0.00001;
while(abs(M - cubic_spline(x, year_values, all_pop)) > 10)
    x = x + h;
end
target = round(cubic_spline(x, year_values, all_pop));
year = round(x, 2);
format longG
month = round((abs(year - fix(year)) * 12));
formattedStr = sprintf("The population will reach %i in %i, %i-th month (Cubic spline)", M, fix(year), month);
disp(formattedStr);

disp(" ")
% ------Task 6 and 7--------
year = 2026;
pred_2026_rural = cubic_spline(year, year_values, rural_pop);
formattedStr = sprintf("Predicted rural population for year %i, is %.0f (Cubic spline)", year, pred_2026_rural);
format LongG
disp(formattedStr);
pred_2026_urban = cubic_spline(year, year_values, urban_pop);
formattedStr = sprintf("Predicted population for year %i, is %.0f (Cubic spline)", year, pred_2026_urban);
format LongG
disp(formattedStr);