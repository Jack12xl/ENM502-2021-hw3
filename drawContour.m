function[] = drawContour(U, title_str)
    figure();
    contourf(U);
    colorbar;
    title(title_str);
end