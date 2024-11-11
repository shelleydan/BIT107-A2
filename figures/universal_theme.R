#Setting the theme for standardization
theme_set(theme_classic(base_size = 15) +
            theme(
              axis.title.y = element_text(face = 'bold', margin = margin(0,20,0,0), size = rel(1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))