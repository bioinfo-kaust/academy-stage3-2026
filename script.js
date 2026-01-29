function openDay(evt, dayName) {
    var i, tabcontent, tabbuttons;

    tabcontent = document.getElementsByClassName("tab-content");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].classList.add("hidden");
    }

    tabbuttons = document.getElementsByClassName("tab-button");
    for (i = 0; i < tabbuttons.length; i++) {
        tabbuttons[i].classList.remove("active");
    }

    const targetTabContent = document.getElementById(dayName);
    targetTabContent.classList.remove("hidden");
    evt.currentTarget.classList.add("active");

    // Check if the content needs to be loaded dynamically for labs
    if (window.location.pathname.includes('labs.html')) {
        const labFileName = dayName + '.html'; // This should be relative to labs.html
        console.log("Attempting to fetch: ", labFileName); // Debug log
        fetch(labFileName)
            .then(response => {
                console.log("Fetch response status: ", response.status); // Debug log
                if (!response.ok) {
                    throw new Error(`HTTP error! status: ${response.status}`);
                }
                return response.text();
            })
            .then(data => {
                const parser = new DOMParser();
                const doc = parser.parseFromString(data, 'text/html');
                // The content we want is within a section with id matching dayName (e.g., labs-day1)
                const labSection = doc.querySelector('#' + dayName);
                if (labSection) {
                    targetTabContent.innerHTML = labSection.innerHTML;
                } else {
                    targetTabContent.innerHTML = '<p>Content not found in the fetched file.</p>';
                }
            })
            .catch(error => {
                console.error('Error loading lab content:', error);
                targetTabContent.innerHTML = '<p>Error loading content.</p>';
            });
    } else if (window.location.pathname.includes('program-schedule.html')) {
        // For program schedule, content is already in the HTML
        // No dynamic loading needed here.
    }
}

// Open the first day by default
document.addEventListener('DOMContentLoaded', (event) => {
    // Check if on the program-schedule.html or labs.html page
    const currentPath = window.location.pathname;

    if (currentPath.includes('program-schedule.html') || currentPath.includes('labs.html')) {
        const firstTabButton = document.querySelector('.tab-buttons .tab-button');
        if (firstTabButton) {
            firstTabButton.click();
        }
    }
});