% for mice whose names are numbers, e.g., '609', so the fieldame is valid
function mouse_field = get_mouse_field(mouse)
    mouse_field = mouse;
    if isstrprop(mouse(1),'digit')
        mouse_field = ['m' mouse_field];
    end
end